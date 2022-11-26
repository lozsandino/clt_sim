
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
            min = 50,
            max = 1e+6,
            step = 50
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
          t(
            paste(
              "The histograms show the distribution of the simulated samples.",
              "The density corresponds to that of the distribution used to",
              "simulate the samples."
            )
          ),
          uiOutput("estVars"),
          plotOutput("hist")
        ),
        tabPanel(
          "Convergence", 
          t(
            paste(
              "Random subsamples of increasing sizes are drawn from the", 
              "simulated samples. The mean difference and its variance",
              "are calculated on each subsample. These plots ilustrate their",
              "convergence and the real values for comparison."
            )
          ),
          uiOutput("estDiff"),
          plotOutput("conv")
        ),
        tabPanel(
          "Bootstrap",
          t(
            paste(
              "This section shows the distribution of the bootstrap mean",
              "difference (out of a 1,000 samples) and compares it to the",
              "normal distribution. Additionally, the result of a Shapiro-Wilk",
              "test (null hypothesis: the bootstrap sample is normally",
              "distributed). Regardless of the result of the normality test,",
              "it's important to evaluate to what extent do the violation",
              "of the assumption may compromise the mean difference parametric",
              "test. The histogram and QQ plot may assist with the latter."
            )
          ),
          uiOutput("shapiroW"),
          plotOutput("boots")
        ),
        tabPanel(
          "Minimum Detectable Difference",
          t(
            paste(
              "Under the assumption of normality, the following procedure can",
              "used to determine the difference that could be detected with a",
              "given sample size, and with some specified significance level",
              "and statistical power (in this case 5% and 80%, respectively).",
              "The orange density corresponds to the distribution of the",
              "difference under the null hypothesis."
            )
          ),
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
      withProgress(
        message = "Processing",
        value = 0,
        {
          # ---- Simulate x and y ----
          setProgress(
            value = 0.1,
            detail = "simulating samples"
          )
          
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
          setProgress(
            value = 0.2,
            detail = "running bootstrap"
          )
          
          boots <- bootsAvgDiff(x, y, 1000)
          densBoots <- dnorm(dSeq(boots), mean(boots), sd(boots))
          swTest <- shapiro.test(boots)
          
          # minimum detectable difference
          mdd <- getMDD(sqrt(realDiffVar), n)
          
          # ---- Plots ----
          ## histograms:
          setProgress(
            value = 0.8,
            detail = "creating plots"
          )
          
          output$hist <- renderPlot({
            p <- grid.arrange(
              ggHist(x, meanX, densX, "x"),
              ggHist(y, meanY, densY, "y"),
              nrow = 1
            )
            
            p
          })
          
          output$estVars <- renderUI({
            withMathJax(
              helpText(
                paste(
                  "$$E(X) =", round(meanX, digits = 2), "; \\quad",
                  "Var(X) =", round(varX, digits = 2), "; \\quad",
                  "\\bar{X} =", round(mean(x), digits = 2), "; \\quad",
                  "{s^2}_X =", round(var(x), digits = 2), "$$"
                )
              ),
              helpText(
                paste(
                  "$$E(Y) =", round(meanY, digits = 2), "; \\quad",
                  "Var(Y) =", round(varY, digits = 2), "; \\quad",
                  "\\bar{Y} =", round(mean(y), digits = 2), "; \\quad",
                  "{s^2}_Y =", round(var(y), digits = 2), "$$"
                )
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
          
          output$estDiff <- renderUI({
            withMathJax(
              helpText(
                paste(
                  "$$E(X) - E(Y) =", round(realDiff, digits = 2), "$$"
                )
              ),
              helpText(
                paste(
                  "$$Var(X - Y) =", round(realDiffVar, digits = 2), "$$"
                )
              )
            )
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
          
          output$shapiroW <- renderUI({
            withMathJax(
              helpText(
                paste(
                  "$$ \\mathrm{Shapiro\\mbox{-}Wilk} \\",
                  "W = ", round(swTest$statistic, digits = 2), "$$"
                )
              ),
              helpText(
                paste(
                  "$$ \\mathrm{p\\mbox{-}value} = ", round(swTest$p.value, digits = 2), "$$"
                )
              )
            )
          })
          
          # minimum detectable difference
          output$plotMDD <- renderPlot({
            ggMDD(x, y, boots)
          })
          
          output$estMDD <- renderUI({
            withMathJax(
              helpText(
                paste(
                  "$$\\mathrm{min. \\ detectable \\ diff.} =", round(mdd, digits = 2), "$$"
                )
              ),
              helpText(
                paste(
                  "$$ n = n_X + n_Y =", 2 * n, "$$"
                )
              )
            )
          })
        }
      )
      
    }
  )
  
}
