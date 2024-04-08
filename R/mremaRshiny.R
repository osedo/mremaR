#'
#' Rshiny app to visualize the mrema() results
#'
#' @export

mremApp <- function(){

  server <- function(input, output, session) {

    results <- shiny::reactive({
      inFile <- input$file
      if (is.null(inFile)) {
        d <- NULL
      } else {
        d <- readRDS(inFile$datapath)
      }
      d
    })

    #print(names(results()))
    shiny::observe({
      shiny::updateSelectInput(session, "CellType",
                        label = "CellType",
                        choices = names(results()),
                        selected = names(results())[1])
    })
    shiny::observe({
      shiny::updateSliderInput(session, "numrow",
                        label = "numrow",
                        max = length(results()))
    })

    table <- shiny::reactive(results()[[input$CellType]]$results %>% DT::datatable(selection = 'single') %>% DT::formatSignif(columns = c("Estimated.Difference", "PVAL", "ADJ.PVAL"), digits = 4))
    output$dTable <-  DT::renderDT(table())
    all.ct <- shiny::reactive(input$across)
    export <- shiny::reactive(input$export)
    central <- shiny::reactive(input$central)
    numrow <- shiny::reactive(input$numrow)
    mixture <- shiny::reactive({
      shiny::req(input$dTable_rows_selected, input$CellType)
      if(all.ct()){
        gs <- results()[[input$CellType]]$results$GENE.SET[input$dTable_rows_selected]
        x <- seq(input$Xrange[1], input$Xrange[2], by = 0.01)
        y <- lapply(seq_along(results()), function(ct) {
          set <- as.numeric(c(results()[[ct]]$parameters[rownames(results()[[ct]]$parameters) == gs,]))
          if(central()){
            y_geneset <- stats::dnorm(x, set[2], sqrt(set[5]))*set[8] + stats::dnorm(x, set[3], sqrt(set[6]))*set[9]
          } else {
            y_geneset <- stats::dnorm(x, set[1], sqrt(set[4]))*set[7] + stats::dnorm(x, set[2], sqrt(set[5]))*set[8] + stats::dnorm(x, set[3], sqrt(set[6]))*set[9]
          }
          data.frame("LFC" = c(x), "Density" = c(y_geneset), "Set" = c(rep(names(results())[ct], length(x))))
        })
        y <- do.call("rbind", y)
      } else {
        gs <- results()[[input$CellType]]$results$GENE.SET[input$dTable_rows_selected]
        x <- seq(input$Xrange[1], input$Xrange[2], by = 0.01)
        set <- as.numeric(c(results()[[input$CellType]]$parameters[rownames(results()[[input$CellType]]$parameters) == gs,]))

        if(central()){
          y_geneset <- stats::dnorm(x, set[2], sqrt(set[5]))*set[8] + stats::dnorm(x, set[3], sqrt(set[6]))*set[9]
          y_background <- stats::dnorm(x, set[11], sqrt(set[14]))*set[17] + stats::dnorm(x, set[12], sqrt(set[15]))*set[18]
        } else {
          y_geneset <- stats::dnorm(x, set[1], sqrt(set[4]))*set[7] + stats::dnorm(x, set[2], sqrt(set[5]))*set[8] + stats::dnorm(x, set[3], sqrt(set[6]))*set[9]
          y_background <- stats::dnorm(x, set[10], sqrt(set[13]))*set[16] + stats::dnorm(x, set[11], sqrt(set[14]))*set[17] + stats::dnorm(x, set[12], sqrt(set[15]))*set[18]
        }
        y <- data.frame("LFC" = c(x, x), "Density" = c(y_geneset, y_background), "Set" = c(rep(gs, length(x)), rep("Background genes", length(x))))

      }
      y
    })
    dist.plot <- shiny::reactive({

      if(numrow() != 0){
        p <- ggplot2::ggplot(data = mixture(), ggplot2::aes(LFC, Density, group = Set)) +
          # geom_line() +
          ggplot2::xlab("Log Fold Change") +
          ggplot2::ylab("Density") +
          ggplot2::geom_area(ggplot2::aes( group = Set, fill = Set), alpha = 0.6, position = "identity") +
          #scale_color_manual(values = natparks.pals(name = "SmokyMtns", n = length(unique(mixture()$Set)), type = "discrete")) +
          ggplot2::scale_fill_brewer(palette = input$pal2) +
          ggplot2::theme_minimal() +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                legend.position = "bottom") +
          ggplot2::facet_wrap(ggplot2::vars(Set), nrow = numrow())
      } else {
        p <- ggplot2::ggplot(data = mixture(), ggplot2::aes(LFC, Density, group = Set)) +
          # geom_line() +
          ggplot2::xlab("Log Fold Change") +
          ggplot2::ylab("Density") +
          ggplot2::geom_area(ggplot2::aes( group = Set, fill = Set), alpha = 0.6, position = "identity") +
          #scale_color_manual(values = natparks.pals(name = "SmokyMtns", n = length(unique(mixture()$Set)), type = "discrete")) +
          ggplot2::scale_fill_brewer(palette = input$pal2) +
          ggplot2::theme_minimal() +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                         legend.position = "bottom")
      }
      p
    })

    output$distPlot <- plotly::renderPlotly({
      dist.plot()
    })

    output$export = shiny::downloadHandler(
      filename = function() {"Plot output.pdf"},
      content = function(file) {
        grDevices::pdf(file)
        plot(dist.plot())
        grDevices::dev.off()
      }
    )

  }

  ui <-  shiny::fluidPage(
     shiny::fluidRow(
       shiny::column(2,  shiny::fluidRow( shiny::fileInput("file", "File input"),
                          shiny::selectInput("CellType", "Cell Type", "")),
             shiny::checkboxInput("across", "Show across cell types", FALSE),
             shiny::checkboxInput("central", "Hide central component", FALSE),
             shiny::sliderInput("Xrange", "X-range:", min = -10, max = 10, value = c(-5, 5)),
             shiny::sliderInput("numrow",label =  "Number of rows:", min = 0, max = 10, value = 0),
             esquisse::palettePicker(
               inputId = "pal2",
               label = "Choose a palette:",
               choices = list(
                 "Set2" = scales::brewer_pal(palette = "Set2")(8),
                 "Paired" = scales::brewer_pal(palette = "Paired")(8),
                 "Dark2" = scales::brewer_pal(palette = "Dark2")(8)
               )
             ),
             shiny::downloadButton('export', "Download pdf")
      ),
       shiny::column(10,  shiny::fluidRow(DT::DTOutput("dTable")),
             plotly::plotlyOutput("distPlot"))
    )
  )


  shiny::shinyApp(ui = ui, server = server)
}

utils::globalVariables(c("LFC", "Density", "Set", "simulation.parameters"))
