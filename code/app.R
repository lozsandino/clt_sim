# rm(list = ls())
# 
# library(shiny)
# library(ggplot2)
# library(gridExtra)
# library(shinybusy)
# 
# source("code/aux_functions.R")

ui <- fluidPage(
  withMathJax(),
  theme = bslib::bs_theme(
    bg = "#e0e0e0",
    fg = "#3b3b3b",
    base_font = "Computer Modern",
    primary = "#3b3b3b"
  ),
  
  add_loading_state(
    ".shiny-plot-output",
    text = "Please wait...",
    svgColor = "steelblue"
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      fluidRow(
        # Distributions
        column(
          6,
          t("Distribution:"),
          selectInput(
            "dst",
            label = "",
            choices = c(
              "normal",
              "gamma"
            )
          )
          
        ),
        
        # Simulation Parameters
        column(
          6,
          
          actionButton(
            "sim",
            label = "Simulate"
          ), 
          align = "center",
          # style = "margin-bottom: 10px;",
          style = "margin-top: 50px;"
        )
      ),
      
      # Parameters
      t("Parameters:"),
      fluidRow(
        column(
          6,
          numericInput(
            "n",
            label = "$$n_X = n_Y$$",
            value = 1000,
            min = 30
          )
        ),
        # Expected value of each distribution
        column(
          6,
          conditionalPanel(
            condition = "input.dst == 'normal'",
            t("$$E(X) = \\mu$$"),
            t("$$Var(X) = \\sigma^2$$")
          ),
          conditionalPanel(
            condition = "input.dst == 'gamma'",
            t("$$E(X) = \\frac{\\alpha}{\\beta}$$"),
            t("$$Var(X) = \\frac{\\alpha}{\\beta^2}$$")
          ),
          style = "margin-top: 15px;"
        )
      ),
      
      ## Normal
      conditionalPanel(
        condition = "input.dst == 'normal'",
        fluidRow(
          column(
            6,
            numericInput(
              "muX",
              label = "$$\\mu_X$$",
              value = 0
            ),
            numericInput(
              "sigmaX",
              label = "$$\\sigma_X$$",
              value = 1,
              min = 0.00001
            )
          ),
          column(
            6,
            numericInput(
              "muY",
              label = "$$\\mu_Y$$",
              value = 0
            ),
            numericInput(
              "sigmaY",
              label = "$$\\sigma_Y$$",
              value = 1,
              min = 0.00001
            )
          )
        )
      ),
      
      ## Gamma
      conditionalPanel(
        condition = "input.dst == 'gamma'",
        fluidRow(
          column(
            6,
            numericInput(
              "alphaX",
              label = "$$\\alpha_x$$",
              value = 1,
              min = 0.00001
            ),
            numericInput(
              "betaX",
              label = "$$\\beta_x$$",
              value = 1,
              min = 0.00001
            )
          ),
          column(
            6,
            numericInput(
              "alphaY",
              label = "$$\\alpha_y$$",
              value = 1,
              min = 0.00001
            ),
            numericInput(
              "betaY",
              label = "$$\\beta_y$$",
              value = 1,
              min = 0.00001
            )
          )
        )
      )
    ),
    
    # Output Plots
    mainPanel(
      width = 9,
      tabsetPanel(
        br(),
        tabPanel(
          "Simulated Samples",
          uiOutput("estX"),
          uiOutput("estY"),
          plotOutput("hist")
        ),
        tabPanel(
          "Convergence", 
          uiOutput("estRealDiff"),
          uiOutput("estSampleDiff"),
          plotOutput("conv")
        ),
        tabPanel(
          "Bootstrap", 
          plotOutput("boots")
        ),
        tabPanel(
          "Minimum Detectable Difference",
          uiOutput("estMDD"),
          plotOutput("plotMDD")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  observeEvent(
    max(input$sim, 1),
    {
      # ---- Simulate x and y ----
      n <- input$n
      
      if (input$dst == "normal") {
        set.seed(NULL)
        
        # x:
        muX <- input$muX
        sigmaX <- input$sigmaX
        
        x <- rnorm(n, muX, sigmaX)
        meanX <- muX
        varX <- sigmaX^2
        
        densX <- dnorm(dSeq(x), muX, sigmaX)
        
        # y:
        muY <- input$muY
        sigmaY <- input$sigmaY
        
        y <- rnorm(n, muY, sigmaY)
        meanY <- muY
        varY <- sigmaY^2
        
        densY <- dnorm(dSeq(y), muY, sigmaY)
        
      } else if (input$dst == "gamma") {
        set.seed(NULL)
        
        # x:
        alphaX <- input$alphaX
        betaX <- input$betaX
        
        x <- rgamma(n, alphaX, betaX)
        meanX <- alphaX/betaX
        varX <- alphaX/betaX^2
        
        densX <- dgamma(dSeq(x), alphaX, betaX)
        
        # y:
        alphaY <- input$alphaY
        betaY <- input$betaY
        
        y <- rgamma(n, alphaY, betaY)
        meanY <- alphaY/betaY
        varY <- alphaY/betaY^2
        
        densY <- dgamma(dSeq(y), alphaY, betaY)
        
      }
      
      # real parameters
      realDiff <- meanX - meanY
      realDiffVar <- varX + varY # Due to independence.
      
      # Bootstrap average difference:
      boots <- bootsAvgDiff(x, y, 1000)
      densBoots <- dnorm(dSeq(boots), mean(boots), sd(boots))
      
      # minimum detectable difference
      mdd <- getMDD(sqrt(realDiffVar), n)
      
      # ---- Plots ----
      ## histograms:
      output$hist <- renderPlot({
        p <- grid.arrange(
          ggHist(x, meanX, densX, "x"),
          ggHist(y, meanY, densY, "y"),
          nrow = 1
        )
        
        p
      })
      
      output$estX <- renderUI({
        withMathJax(
          paste(
            "$$E(X) =", round(meanX, digits = 2), "; \\quad",
            "Var(X) =", round(varX, digits = 2), "; \\quad",
            "\\bar{X} =", round(mean(x), digits = 2), "; \\quad",
            "{s^2}_X =", round(var(x), digits = 2), "$$"
          )
        )
      })
      
      output$estY <- renderUI({
        withMathJax(
          paste(
            "$$E(Y) =", round(meanY, digits = 2), "; \\quad",
            "Var(Y) =", round(varY, digits = 2), "; \\quad",
            "\\bar{Y} =", round(mean(y), digits = 2), "; \\quad",
            "{s^2}_Y =", round(var(y), digits = 2), "$$"
          )
        )
      })
      
      ## convergence of parameters:
      output$conv <- renderPlot({
        p <- grid.arrange(
          ggConv(x - y, realDiff),
          ggConv(x - y, realDiffVar, var = TRUE),
          nrow = 1
        )
        
        p
      })
      
      ## bootstrap:
      output$boots <- renderPlot({
        p <- grid.arrange(
          ggHist(boots, realDiff, densBoots, nameX = "Boostrap Avg. Diff."),
          ggNorm(boots),
          widths = c(2, 1),
          nrow = 1
        )
        
        p
      })
      
      output$estRealDiff <- renderUI({
        withMathJax(
          paste(
            "$$E(X) - E(Y) =", round(realDiff, digits = 2), "; \\quad", 
            "Var(X - Y) =", round(realDiffVar, digits = 2), "$$"
          )
        )
      })
      
      output$estSampleDiff <- renderUI({
        withMathJax(
          paste(
            "$$\\bar{X} - \\bar{Y} =", round(mean(x) - mean(y), digits = 2),"; \\quad", 
            "{s^2}_{X-Y} =", round(var(x) + var(y), digits = 2), "$$"
          )
        )
      })
      
      # minimum detectable difference
      output$plotMDD <- renderPlot({
        ggMDD(x, y, boots)
      })
      
      output$estMDD <- renderUI({
        withMathJax(
          paste(
            "$$mdd =", round(mdd, digits = 2), "$$"
          )
        )
      })
      
    }
  )
  
}

# shinyApp(ui, server,  options = list(width = 1200, height = 580))
