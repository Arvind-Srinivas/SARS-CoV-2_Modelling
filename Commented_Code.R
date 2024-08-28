#Load necessary libraries.
library(deSolve)
library(ggplot2)
library(dplyr)
library(viridis)
library(shiny)


#---- DEFINING SIR-ODES AND CONDITIONS ----

#Defining the ODE system. 4 compartmental SIR model with two infectious compartments symbolizing short-term and persistent strains.
sirodes <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta_s <- R0s * gamma_s
    beta_p <- R0p * gamma_p
    dS <- (mu * N) - (beta_s * S * I_s / N) - (beta_p * S * I_p / N) - (mu * S) + (rho * R)
    dI_s <- (beta_s * S * I_s / N) - (gamma_s * I_s) - (mu * I_s)
    dI_p <- (beta_p * S * I_p / N) - (gamma_p * I_p) - (mu * I_p)
    dR <- (gamma_s * I_s) + (gamma_p * I_p) - (rho * R) - (mu * R)
    list(c(dS, dI_s, dI_p, dR))
  })
}

#Constant parameters. Refer to table 1 in the written report.
mu <- 2.5e-5
N <- 2794356
gamma_s <- 0.048
gamma_p <- 0.0056

#Initial conditions. Population of Toronto was used as the initial population number. 1000 Is infections and 1 persistent infection at the start of the simulation. 
initial_state <- c(S = N - 1001, I_s = 1000, I_p = 1, R = 0)

#Time points for simulating (1825 days or 5 years).
times <- seq(0, 1825, by = 1)

#---- PLOT: Final Number of Persistent Infections vs. R0p ----

#This population dynamics UI was created to observe how the model performs under realistic and unrealistic conditions for R0p. The model performs well for outlier values of R0p but not for values within the R0p range. With higher R0p values, the number of persistent infections exponentially increase, which contradicts the gradual increase of persistent infections in real-time. 

#Define a range for R0p. A large range was chosen to show the limitations of the model. R0p must be between 0-4 as the assumption for persistent infections is low transmissibility. 
R0p_values <- seq(0.1, 50, by = 0.5)

#Initialize a vector to store the final number of persistent infections.
final_Ip_values <- numeric(length(R0p_values))

#Run the model for each R0p value. R0s and immunity waning (p) were constants with values within their acceptable ranges. 
for (i in 1:length(R0p_values)) {
  parameters <- c(
    R0s = 8,
    R0p = R0p_values[i],
    gamma_s = gamma_s,
    gamma_p = gamma_p,
    mu = mu,
    N = N,
    rho = 0.005
  )
  
  #Running the model.
  out <- ode(y = initial_state, times = times, func = sirodes, parms = parameters)
  out <- as.data.frame(out)
  
  #Recording the final number of persistent infections.
  final_Ip_values[i] <- tail(out$I_p, 1)
}

#Create a data frame for plotting.
df_final_Ip <- data.frame(R0p = R0p_values, Final_Ip = final_Ip_values)

#Data exploration of dataframe: 
head(df_final_Ip)
summary(df_final_Ip)
sum(is.na(df_final_Ip))

#Plotting the final number of persistent infections vs. R0p.
plot_final_Ip <- ggplot(df_final_Ip, aes(x = R0p, y = Final_Ip)) +
  geom_line(size = 1) + 
  geom_vline(xintercept = 3, linetype = "dashed", color = "orange", size = 0.8) +  
  labs(title = "Final Number of Persistent Infections vs. R0p",
       x = "Persistent Reproduction Value (R0p)",
       y = "Final Number of Persistent Infections (Ip)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(plot_final_Ip)


#---- UI INTEGRATION FOR POPULATION DYNAMICS ----

#Define UI for the app. There are 3 sliders for variable parameters: short-term reproduction number, persistent reproduction number, and immune waning rate. 
ui <- fluidPage(
  titlePanel("SIR Model: Population Dynamics"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("R0s", "R0s (Short-term Reproduction Number)", min = 1, max = 20, value = 2, step = 0.5),
      sliderInput("R0p", "R0p (Persistent Reproduction Number)", min = 0.1, max = 20, value = 0.5, step = 0.1),
      sliderInput("rho", "Immunity Waning (rho)", min = 0.001, max = 0.01, value = 0.005, step = 0.001)
    ),
    
    mainPanel(
      plotOutput("populationPlot")
    )
  )
)

#Defining server logic.
server <- function(input, output) {
  
  #Reactive expression to update the ODE system based on UI inputs.
  populationData <- reactive({
    
    #Parameters for the ODE.
    parameters <- c(
      R0s = input$R0s,
      R0p = input$R0p,
      gamma_s = gamma_s,
      gamma_p = gamma_p,
      mu = mu,
      N = N,
      rho = input$rho
    )
    
    #Run the ODE system.
    out <- ode(y = initial_state, times = times, func = sirodes, parms = parameters)
    as.data.frame(out)
  })
  
  #Plot the population dynamics based on UI inputs.
  output$populationPlot <- renderPlot({
    out <- populationData()
    
    #Colorblind-friendly palette.
    cb_palette <- c("Susceptible" = "#0072B2", "Short-term Infections" = "#D55E00", 
                    "Persistent Infections" = "#009E73", "Recovered" = "#CC79A7")
    
    #Plotting.
    ggplot(out, aes(x = time)) +
      geom_line(aes(y = S, color = "Susceptible"), size = 1.2) +
      geom_line(aes(y = I_s, color = "Short-term Infections"), size = 1.2) +
      geom_line(aes(y = I_p, color = "Persistent Infections"), size = 1.2) +
      geom_line(aes(y = R, color = "Recovered"), size = 1.2) +
      labs(title = "Population Dynamics",
           x = "Time (days)",
           y = "Number of Infections",
           color = "Compartment") +
      scale_color_manual(values = cb_palette) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  })
}

#Run the application.
shinyApp(ui = ui, server = server)

#The two graphs in the written report represent specific R0p values that are realistic and unrealistic. Adjusting these sliders for any R0p value would show how values much greater than the 0-4 range would show an exponential increase of persistent infections. To replicate the graphs: 

#Realistic Conditions: set R0p = 2.5, R0s = 8, immunity waning (p) = 0.005

#Unrealistic Conditions: set R0p = 20, R0s = 8, immunity waning (p) = 0.005


#---- SENSITIVITY ANALYSIS: rho (p) vs R0p ----

#Define ranges for R0p and rho. I wanted to show how unrealistic values of R0p influence the number of persistent infections, so I went beyond the 0-4 range for R0p_values.
R0p_values <- seq(0.1, 20, by = 0.5)  
rho_values <- seq(0.001, 0.01, by = 0.001)

#Create a data frame to store results.
results_rho_R0p <- expand.grid(R0p = R0p_values, rho = rho_values, final_Ip = NA)

#Perform the multi-sensitivity analysis for rho vs R0p. R0s was constant at 8. 
for (i in 1:nrow(results_rho_R0p)) {
  parameters <- c(
    R0s = 8,
    R0p = results_rho_R0p$R0p[i],
    gamma_s = gamma_s,
    gamma_p = gamma_p,
    mu = mu,
    N = N,
    rho = results_rho_R0p$rho[i]
  )
  
  #Run the model.
  out <- ode(y = initial_state, times = times, func = sirodes, parms = parameters)
  out <- as.data.frame(out)
  
  #Record the final number of persistent infections.
  results_rho_R0p$final_Ip[i] <- tail(out$I_p, 1)
}

#Visualizing the results using a heatmap for rho vs R0p.
heatmap_plot_rho_R0p <- ggplot(results_rho_R0p, aes(x = R0p, y = rho, fill = final_Ip)) +
  geom_tile() +
  scale_fill_viridis(option = "plasma") +  
  labs(title = "Final Number of Persistent Infections (Ip) Sensitivity Analysis: Immunity Waning (rho) vs R0p",
       x = "Persistent Reproduction Number (R0p)",
       y = "Immunity Waning (rho)",
       fill = "Final Ip") +
  theme_minimal()

print(heatmap_plot_rho_R0p) 


#---- SENSITIVITY ANALYSIS: R0s vs R0p ----

#Define ranges for R0s and R0p. I chose a large range to show how outliers allow the model to perform as expected.
R0s_values <- seq(1, 20, by = 0.5)  
R0p_values <- seq(0.1, 20, by = 0.5)

#Create a data frame to store results.
results_R0s_R0p <- expand.grid(R0s = R0s_values, R0p = R0p_values, final_Ip = NA)

#Perform the multisensitivity analysis for R0s vs R0p. Rho is constant at a value of 0.005.
for (i in 1:nrow(results_R0s_R0p)) {
  parameters <- c(
    R0s = results_R0s_R0p$R0s[i],
    R0p = results_R0s_R0p$R0p[i],
    gamma_s = gamma_s,
    gamma_p = gamma_p,
    mu = mu,
    N = N,
    rho = 0.005  
  )
  
  #Run the model.
  out <- ode(y = initial_state, times = times, func = sirodes, parms = parameters)
  out <- as.data.frame(out)
  
  #Record the final number of persistent infections.
  results_R0s_R0p$final_Ip[i] <- tail(out$I_p, 1)
}

#Visualizing the results using a heatmap for R0s vs R0p.
heatmap_plot_R0s_R0p <- ggplot(results_R0s_R0p, aes(x = R0s, y = R0p, fill = final_Ip)) +
  geom_tile() +
  scale_fill_viridis(option = "plasma") +  
  labs(title = "Final Number of Persistent Infections (Ip) Sensitivity Analysis: R0s vs R0p",
       x = "Short-term Reproduction Number (R0s)",
       y = "Persistent Reproduction Number (R0p)",
       fill = "Final Ip") +
  theme_minimal()

print(heatmap_plot_R0s_R0p)


#---- SENSITIVITY ANALYSIS: R0s vs rho ----

#Define ranges for R0s and rho.
R0s_values <- seq(1, 20, by = 0.5)  
rho_values <- seq(0.001, 0.01, by = 0.001)

#Create a data frame to store results.
results_R0s_rho <- expand.grid(R0s = R0s_values, rho = rho_values, final_Ip = NA)

#Perform the multisensitivity analysis for R0s vs rho. R0p is constant at 2.5. 
for (i in 1:nrow(results_R0s_rho)) {
  parameters <- c(
    R0s = results_R0s_rho$R0s[i],
    R0p = 2.5,  
    gamma_s = gamma_s,
    gamma_p = gamma_p,
    mu = mu,
    N = N,
    rho = results_R0s_rho$rho[i]
  )
  
  #Run the model.
  out <- ode(y = initial_state, times = times, func = sirodes, parms = parameters)
  out <- as.data.frame(out)
  
  #Record the final number of persistent infections.
  results_R0s_rho$final_Ip[i] <- tail(out$I_p, 1)
}

#Visualizing the results using a heatmap for R0s vs rho.
heatmap_plot_R0s_rho <- ggplot(results_R0s_rho, aes(x = R0s, y = rho, fill = final_Ip)) +
  geom_tile() +
  scale_fill_viridis(option = "plasma") +  #Use plasma color scale.
  labs(title = "Final Number of Persistent Infections (Ip) Sensitivity Analysis: R0s vs Immunity Waning (rho)",
       x = "Short-term Reproduction Number (R0s)",
       y = "Immunity Waning (rho)",
       fill = "Final Ip") +
  theme_minimal()

print(heatmap_plot_R0s_rho)
