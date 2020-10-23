library(shiny)
library(tidyverse)
library(igraph)
library(RColorBrewer)
library(ggpubr)

ui <- fluidPage(
  titlePanel("Disease modeling and graph theory"),
  sidebarLayout(
    # Sidebar panel for inputs
    sidebarPanel(
      numericInput(
        inputId = "replicates",
        label = "Number of replicates (10-100)",
        value = 20,
        min = 10,
        max = 100,
        step = 5
        ), # close numericInput
      numericInput(
        inputId = "n",
        label = "Number of individuals (20-100)",
        value = 50,
        min = 20,
        max = 100,
        step = 5
        ),
      numericInput(
        inputId = "days",
        label = "Length of simulation (10-100)",
        value = 20,
        min = 10,
        max = 100,
        step = 5
      ),
      h4("Beta distribution parameters for edge weights"),
      sliderInput(
        inputId = "alpha",
        label = "Alpha",
        min = 1,
        max = 10,
        step = 1,
        value = 1
      ),
      sliderInput(
        inputId = "beta",
        label = "Beta",
        min = 1,
        max = 10,
        step = 1,
        value = 9
      ),
      h4("Probability of presence of edge"),
      sliderInput(
        inputId = "edgeprob",
        label = "Probability",
        min = 0.01,
        max = 0.99,
        step = 0.01,
        value = 0.1
      ),
      actionButton("go", label = "Run Simulation")
      ), # close sidebarPanel
    mainPanel(
      textOutput("introduction"),
      uiOutput("github_url"),
      h3("Infection summary"),
      plotOutput(outputId = "sim1"),
      h3("Edge weights"),
      textOutput(outputId = "part2"),
      plotOutput(outputId = "sim2"),
      h3("Graph"),
      textOutput(outputId = "part3"),
      plotOutput(outputId = "igraph")

    )
    ) # close sidebarLayout

  
) # close fluidPage

server <- function(input, output) {
  
  ## URL for my github page
  github_url <- a("Click here for code on my github page", href = "https://github.com/mgdesaix/disease-modeling")
  output$github_url <- renderUI({
    tagList(github_url)
  })
  
  colors2 <- brewer.pal(3, "Dark2")
  
  disease_sim <- function(n, replicates, day, edge.prob, sim.type, alpha, beta){
    n.edges <- sum(1:(n-1))
    
    final.daily.tracker <- matrix(nrow = n,
                                  ncol = replicates,
                                  0)
    
    # matrix to be returned of # of infected individuals daily
    # note, rows = days, columns = replicates
    final.daily.infected <- matrix(nrow = day,
                                   ncol = replicates,
                                   0)
    # matrix for R0
    final.infected.number <- matrix(nrow = n,
                                    ncol = replicates,
                                    0)
    # list of all starting adj.mats
    starting.adj.mat.list <- list()
    
    starting.individual <- c()
    ########################### FUNCTIONs####################################
    # function for converting days to appropriate infection weight
    # use with sapply
    get.weights <- function(x){
      if(x == 0){return(0)}
      else if(x >= 1 & x < 4){return(1)}
      else if(x >= 4 & x < 7){return(0.1)}
      else if(x >= 7){return(0)}
    }
    
    change.adj.matrix <- function(input.adj.mat, daily.tracker){
      for(i in 1:length(daily.tracker)){
        if(daily.tracker[i,1] == 1){
          input.adj.mat[,i] <- 0
        }
      }
      return(input.adj.mat)
    }
    
    # make sure an individual infects no more than 1 person per day
    check.infected <- function(infected.mat){
      for(i in 1:ncol(infected.mat)){
        if(sum(infected.mat[,i]) > 1){
          # this is a way to deal with ind infected by multiple others
          infected.index <- sample(which(infected.mat[,i]),1)
          infected.mat[,i] <- FALSE
          infected.mat[infected.index, i] <- TRUE
        }
      }
      return(infected.mat)
    }
    
    # add 1 for each individual infected
    update.infection.tracker <- function(infectors, infection.tracker){
      for(i in 1:length(infectors)){
        infection.tracker[infectors[i],1] <- infection.tracker[infectors[i],1] + 1
      }
      return(infection.tracker)
    }
    
    # add 1 for new infections and old, stop at day 11
    update.daily.tracker <- function(infected.mat, daily.tracker){
      daily.tracker <- sapply(daily.tracker, function(x) 
        if(x > 0 & x < 11){return(x+1)} 
        else {return(x)}) %>%
        matrix(ncol = 1)
      infected.ind <- which(t(infected)) %% n
      for(i in infected.ind){
        daily.tracker[i,1] <- 1
      }
      return(daily.tracker)
    }
    
    # making sure starting individual can infect others
    get.start.ind <- function(adj.mat){
      # sum by rows
      val <- apply(adj.mat, 1, sum)
      # which rows are greater than 0
      start.ind <- which(val > 0) %>%
        sample(size = 1)
    }
    
    ######################################################
    
    for(i in 1:replicates){
      # create a random sample of edge values
      # edge.values.weighted <- runif(n = n.edges, 
      #                               min = 0, 
      #                               max = 1) * rbinom(n = n.edges, 
      #                                                 size = 1, 
      #                                                 prob = edge.prob)
      
      # or change things up with the beta distribution
      edge.values.weighted <- rbeta(n = n.edges, 
                                    shape1 = alpha, 
                                    shape2 = beta) * rbinom(n = n.edges, 
                                                         size = 1, 
                                                         prob = edge.prob)
      
      adj.mat <- matrix(nrow = n, ncol = n, 0)
      # create just a lower triangle
      adj.mat[lower.tri(adj.mat)] <- edge.values.weighted
      # full matrix needed for keeping track of the SUPER SPREADERS
      adj.mat <- adj.mat + t(adj.mat)
      
      ################### Simulation type changes here ####################
      
      if(sim.type == "A"){
        random.50 <- sample(x = c(1:n),
                            size = round(0.5*n),
                            replace = F)
        index <- c(1:n) %in% random.50
        sim.type.A.mat <- matrix(nrow = n,
                                 ncol = 1,
                                 1)
        sim.type.A.mat[index,] <- 0.4
        adj.mat <- adj.mat * sim.type.A.mat[,1]
      } else if(sim.type == "B"){
        random.70 <- sample(x = c(1:n),
                            size = round(0.7 * n),
                            replace = F)
        index <- c(1:n) %in% random.70
        sim.type.B.mat <- matrix(nrow = n,
                                 ncol = 1,
                                 1)
        sim.type.B.mat[index,] <- 0.25
        adj.mat <- adj.mat * sim.type.B.mat[,1]
      } else{
        adj.mat <- adj.mat
      }
      
      
      ###################################################################
      
      # random individual to start the infection
      start.ind <- get.start.ind(adj.mat)
      # create 0 matrix
      daily.tracker <- matrix(nrow = n, ncol = 1, 0)
      # make entire row for individual equal to 1
      daily.tracker[start.ind,] <- 1
      
      # vector of infection tracker
      # i.e. number of individuals infected by each individual
      infection.tracker <- matrix(nrow = n, ncol = 1, 0)
      
      infected.daily <- c()
      
      starting.adj.mat.list[[i]] <- adj.mat
      starting.individual[i] <- start.ind
      
      for(j in 1:day){
        
        infected.daily[j] <- sum(daily.tracker > 0 & daily.tracker < 11)
        # get weights for each individual based on day
        infection.weights <- sapply(daily.tracker, get.weights) %>%
          as.matrix()
        # convert infected individuals probability of infection to 0
        adj.mat <- change.adj.matrix(input.adj.mat = adj.mat,
                                     daily.tracker = daily.tracker)
        
        # roll the dice
        sample.mat <- runif(n = n^2, min = 0, max = 1) %>%
          matrix(nrow = n)
        # who got infected?
        infected <- sample.mat < (adj.mat * infection.weights[,1])
        # reduce to 1 indiv infecting
        infected <- check.infected(infected)
        # who did the infecting?
        infectors <- which(infected) %% n
        # update infection.tracker
        infection.tracker <- update.infection.tracker(
          infectors = infectors,
          infection.tracker = infection.tracker)
        # update daily tracker
        daily.tracker <- update.daily.tracker(infected.mat = infected,
                                              daily.tracker)
      }
      
      final.daily.infected[,i] <- infected.daily
      final.infected.number[,i] <- infection.tracker
      final.daily.tracker[,i] <- daily.tracker
      
    }
    
    out.list <- list("FinalDailyTracker" = final.daily.tracker,
                     "FinalDailyInfected" = final.daily.infected,
                     "FinalInfectedNumber" = final.infected.number,
                     "Replicates" = replicates,
                     "AdjMat" = starting.adj.mat.list[[1]])
    return(out.list)
  }
  
  sim1.input <- eventReactive(input$go, {
    disease_sim(n = input$n,
                replicates = input$replicates,
                day = input$days,
                alpha = input$alpha,
                beta = input$beta,
                edge.prob = input$edgeprob,
                sim.type = "Regular")
  })
  
  sim2.input <- eventReactive(input$go, {
    disease_sim(n = input$n,
                replicates = input$replicates,
                day = input$days,
                alpha = input$alpha,
                beta = input$beta,
                edge.prob = input$edgeprob,
                sim.type = "A")
  })
  
  sim3.input <- eventReactive(input$go, {
    disease_sim(n = input$n,
                replicates = input$replicates,
                day = input$days,
                alpha = input$alpha,
                beta = input$beta,
                edge.prob = input$edgeprob,
                sim.type = "B")
  })
  
  output$sim1 <- renderPlot({
    sim1 <- sim1.input() # input from go button
    sim2 <- sim2.input()
    sim3 <- sim3.input()
    
    sim1.1 <- sim1$FinalDailyInfected %>%
      as_tibble() %>%
      mutate(Daily_mean = rowMeans(.)) %>%
      add_column(Day = 1:input$days) %>%
      select(Day, Daily_mean) %>%
      add_column(Simulation = "Standard")
      
    sim2.1 <- sim2$FinalDailyInfected %>%
      as_tibble() %>%
      mutate(Daily_mean = rowMeans(.)) %>%
      add_column(Day = 1:input$days) %>%
      select(Day, Daily_mean) %>%
      add_column(Simulation = "A")
    
    sim3.1 <- sim3$FinalDailyInfected %>%
      as_tibble() %>%
      mutate(Daily_mean = rowMeans(.)) %>%
      add_column(Day = 1:input$days) %>%
      select(Day, Daily_mean) %>%
      add_column(Simulation = "B")
    
    replicate.r0.1 <- apply(sim1$FinalInfectedNumber, 2, mean) %>%
      as_tibble() %>%
      add_column(group = "Standard")
    replicate.r0.2 <- apply(sim2$FinalInfectedNumber, 2, mean) %>%
      as_tibble() %>%
      add_column(group = "A")
    replicate.r0.3 <- apply(sim3$FinalInfectedNumber, 2, mean) %>%
      as_tibble() %>%
      add_column(group = "B")
    
    
    p2 <- rbind(replicate.r0.1, replicate.r0.2) %>%
      rbind(replicate.r0.3) %>%
      ggplot(aes(x = group, y = value, color = group)) +
      geom_boxplot() +
      theme_bw() +
      scale_color_manual(name = "Model",
                         values = colors2,
                         breaks = c("Standard",
                                    "A",
                                    "B"),
                         labels = c("Standard",
                                    "Health Program A",
                                    "Health Program B")) +
      labs(y = "R0",
           x = "Models") +
      theme(legend.position = "none")
    
    p1 <- rbind(sim1.1, sim2.1) %>%
      rbind(sim3.1) %>%
      ggplot(aes(x = Day, 
                 y = Daily_mean,
                 group = Simulation,
                 color = Simulation
                 )) +
      geom_point() +
      scale_color_manual(name = "Model",
                         values = colors2,
                         breaks = c("Standard",
                                    "A",
                                    "B"),
                         labels = c("Standard",
                                    "Health Program A",
                                    "Health Program B")) +
      labs(x = "Day",
           y = "Number infected daily") +
      theme_bw()
    
    p.combined <- ggarrange(p1, p2, nrow = 1)
    p.combined
  })
  
  output$sim2 <- renderPlot({
    sim1 <- sim1.input() # input from go button
    sim2 <- sim2.input()
    sim3 <- sim3.input()
    colors2 <- brewer.pal(3, "Dark2")
    
    ########## P1 ########################
    
    p1 <- tibble("Standard" = sim1$AdjMat %>% as.vector(),
           "A" = sim2$AdjMat %>% as.vector(),
           "B" = sim3$AdjMat %>% as.vector()) %>%
      pivot_longer(cols = c(Standard, A, B),
                   names_to = "Model") %>% 
      ggplot() + 
      geom_histogram(aes(x = value, fill = Model), 
                     alpha = 0.6, 
                     position = "identity",
                     bins = 20) +
      theme_bw() +
      scale_fill_manual(name = "Model",
                         values = colors2,
                         breaks = c("Standard",
                                    "A",
                                    "B"),
                         labels = c("Standard",
                                    "Health Program A",
                                    "Health Program B")) +
      labs(x = "Edge weights",
           y = "Count")
    
    ########## P2 ########################
    
    p2 <- tibble("Standard" = sim1$AdjMat %>% as.vector(),
                  "A" = sim2$AdjMat %>% as.vector(),
                  "B" = sim3$AdjMat %>% as.vector()) %>%
      pivot_longer(cols = c(Standard, A, B),
                   names_to = "Model") %>% 
      ggplot() + 
      geom_boxplot(aes(x = Model, y = value, fill = Model),
                   aes = 0.6) +
      theme_bw() +
      scale_fill_manual(name = "Model",
                        values = colors2,
                        breaks = c("Standard",
                                   "A",
                                   "B"),
                        labels = c("Standard",
                                   "Health Program A",
                                   "Health Program B")) +
      labs(x = "Model",
           y = "Edge weight")
    
    ########## P3 ########################
    sim1.edges <- sim1$AdjMat %>% apply(., 1, function(x) sum(x > 0))
    sim2.edges <- sim2$AdjMat %>% apply(., 1, function(x) sum(x > 0))
    sim3.edges <- sim3$AdjMat %>% apply(., 1, function(x) sum(x > 0))
  
    
    p3 <- tibble("Standard" = sim1.edges,
                      "A" = sim2.edges,
                      "B" = sim3.edges) %>%
      pivot_longer(cols = c(Standard, A, B),
                   names_to = "Model") %>% 
      ggplot() + 
      geom_boxplot(aes(x = Model, y = value, fill = Model),
                   aes = 0.6) +
      theme_bw() +
      scale_fill_manual(name = "Model",
                        values = colors2,
                        breaks = c("Standard",
                                   "A",
                                   "B"),
                        labels = c("Standard",
                                   "Health Program A",
                                   "Health Program B")) +
      labs(x = "Model",
           y = "Number of edges per individual")
    
    p.combined <- ggarrange(p1, p2, p3, nrow = 1)
    p.combined
  })
  
  output$igraph <- renderPlot({
    sim1 <- sim1.input()
    sim2 <- sim2.input()
    sim3 <- sim3.input()
    
    undir.graph1 <- graph_from_adjacency_matrix(sim1$AdjMat, 
                                               mode = "undir",
                                               weighted = TRUE)
    undir.graph2 <- graph_from_adjacency_matrix(sim2$AdjMat, 
                                                mode = "undir",
                                                weighted = TRUE)
    undir.graph3 <- graph_from_adjacency_matrix(sim3$AdjMat, 
                                                mode = "undir",
                                                weighted = TRUE)
    
    par(mfrow=c(1,3))
    plot(undir.graph1,
         edge.width=(edge_attr(undir.graph1)$weight * 10),
         main = "Standard")
    plot(undir.graph2,
         edge.width=(edge_attr(undir.graph2)$weight * 10),
         main = "A")
    plot(undir.graph3,
         edge.width=(edge_attr(undir.graph3)$weight * 10),
         main = "B")

  })
  
  output$part2 <- renderText(
    paste("Here are overlayed histograms for the starting edge weights in the adjacency matrices for the 3 different models")
  )
  
  output$part3 <- renderText(
    paste("Example of a starting graph with edges thickness proportional to edge weight")
  )
  
  output$introduction <- renderText({
    paste("Here are different disease modeling scenarios based on graph theory.", "
          Events are tracked on a daily basis with a network of individuals. ",
          "Individuals are tracked as nodes with weighted edges",
          "(the weights are the probability of infection – ", 
          "assume that it is the same probability to have node A infect node B as it is B to infect A).", 
          "A random individual is initially infected." ,
          "Each day that individual may infect any adjacent individual with a probability equal to the edge weight for up to 3 days after infection.", 
          "On day 4 the infected individual can infect individuals but now the edge weights are only 10% of the original", 
          "(to reflect that the individual is showing symptoms and mostly isolated).", 
          "After day 6 the individual is no longer able to infect anyone adjacent (consider how to modify the graph appropriately)", 
          "but is still infected. After 10 days the individual is no longer infected (and the immune response is such that the individual", 
          "cannot be reinfected over the period of the simulation.",
          "In one model (Health Program A), 50% of individuals have their edge weights reduced by 60%.",
          "In a second model (Health Program B), 70% of individuals have an edge weight reduction by 80%.",
          "For each model replicate, the average number of individuals every infected individual infected is R0.")
  })
}

shinyApp(ui = ui,
         server = server)