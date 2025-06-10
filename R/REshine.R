#'
#' Rshiny app to visualize the REtest() results
#'
#' @export


REshine <- function(){

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

    table <- shiny::reactive(results()[[input$CellType]]$results[,c(1:3,6)] %>% DT::datatable(selection = 'single') %>% DT::formatSignif(columns = c("Prop.DE", "P.value"), digits = 4))
    output$dTable <-  DT::renderDT(table())
    all.ct <- shiny::reactive(input$across)
    export <- shiny::reactive(input$export)
    #disc <- shiny::reactive(input$disc)
    numrow <- shiny::reactive(input$numrow)
    mixture <- shiny::reactive({
      shiny::req(input$dTable_rows_selected, input$CellType)
      i <- results()[[input$CellType]]$results$Gene.Set[input$dTable_rows_selected]
      x <- seq(input$Xrange[1], input$Xrange[2], by = 0.01) # unique(c(seq(input$Xrange[1], -input$thresh, by = 0.01), seq(input$thresh, input$Xrange[2], by = 0.01)))#
      y <- lapply(seq_along(results()), function(ct) {

        if(i %in% results()[[ct]]$results$Gene.Set){
          i <- which(results()[[ct]]$results$Gene.Set == i)


        inset_index <- results()[[ct]]$results$Index[[i]]
        outset_index <- c(1:nrow(results()[[ct]]$results))[-inset_index]
        y_outset <- unlist(lapply(x, function(x_i){
          mean(stats::dnorm(x_i, results()[[ct]]$estimates$lfc[outset_index], results()[[ct]]$estimates$lfcSE[outset_index]))
        }))
        y_inset <- unlist(lapply(x, function(x_i){
          mean(stats::dnorm(x_i, results()[[ct]]$estimates$lfc[inset_index], results()[[ct]]$estimates$lfcSE[inset_index]))
        }))
        f.max <- max(c(y_inset[which(x == 0) + 1], y_outset[which(x == 0) + 1]))
        s.max <- min(c(max(y_inset), max(y_outset)))
        t.max <- max(c(max(y_inset), max(y_outset)))
        density <- data.frame("LFC" = c(x, x), "Density" = c(y_inset, y_outset), "Set" = c(rep(results()[[ct]]$results$Gene.Set[i], length(x)), rep("Background genes", length(x))), "Cell" = c(rep(names(results())[ct], length(x)*2)), "pval" = rep(results()[[ct]]$results$P.value[i], length(x)*2), "Weight.DE" = rep(results()[[ct]]$results$Prop.DE[i], length(x)*2))
        #density$text <- paste0("P-value: ", signif(density$pval, 4), "<br>", "DE weight: ", signif(density$Weight.DE, 4))#paste0("Set: ", density$Set, "<br>", "Cell Type: ", density$Cell, "<br>", "P-value: ", signif(density$pval, 4), "<br>", "DE weight: ", signif(density$Weight.DE, 4))
        density$Density[which(density$LFC > -input$thresh & density$LFC < input$thresh)] <- 0
        #print(head(density))
        density
        }
      })
      #y <- do.call("rbind", y)

      names(y) <- names(results())
      y
    })



    dist.plot <- shiny::reactive({
      if(all.ct()){
        pd <- do.call("rbind", mixture())
        pl <- ggplot2::ggplot(data = pd, ggplot2::aes(LFC, Density, group = Set, colour = Set)) +
          ggplot2::xlab("Log Fold Change") +
          ggplot2::ylab("Density") +
          ggplot2::geom_area(ggplot2::aes( group = Set, fill = Set), alpha = 0.6, position = "identity") +
          ggplot2::scale_color_brewer(palette = input$pal2) +
          ggplot2::scale_fill_brewer(palette = input$pal2) +
          #ggplot2::annotate("text", x = (max(pd$LFC) - 1), y = max(pd$Density) - max(pd$Density)*0.1, label = paste0("P-value - ", pd$pval)) +
          ggplot2::facet_wrap(ggplot2::vars(Cell), nrow = input$numrow) +
          ggplot2::theme_minimal() +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = "bottom")

        text_data <- data.frame("label" = paste0("P-value - ", unique(signif(pd$pval, 3))), "Cell" = unique(pd$Cell))
        text_data$Set <- pd$Set[1]
        #text_data <- data.frame("label" = unlist(lapply(paste0("P-value -", unique(signif(pd$pval, 4))), function(x) rep(x, nrow(pd)/length(unique(pd$Cell))))), "Cell" = pd$Cell)
        #print(text_data)
        pl <- pl + ggplot2::geom_text(size = 2.5, data = text_data, colour = "black", mapping = ggplot2::aes((max(pd$LFC) - 2), y = max(pd$Density) - max(pd$Density)*0.1, label = label))

      } else {
        pd <- mixture()[[input$CellType]]
        pl <- ggplot2::ggplot(data = pd, ggplot2::aes(LFC, Density, group = Set, colour = Set)) +
          ggplot2::xlab("Log Fold Change") +
          ggplot2::ylab("Density") +
          ggplot2::geom_area(ggplot2::aes( group = Set, fill = Set), alpha = 0.6, position = "identity") +
          ggplot2::scale_color_brewer(palette = input$pal2) +
          ggplot2::scale_fill_brewer(palette = input$pal2) +
          ggplot2::theme_minimal() +
          ggplot2::annotate("text", x = (max(pd$LFC) - 2), y = max(pd$Density) - max(pd$Density)*0.1, label = unique(paste0("P-value - ", signif(pd$pval, 3))), size = 2.5) +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5), legend.position = "bottom")
      }


      pl

    })


    output$distPlot <- shiny::renderPlot({
      dist.plot()
    })

    output$distPlotly <- plotly::renderPlotly({
      plotly::ggplotly(dist.plot()) %>%
        plotly::layout(legend = list(orientation = 'h', x = 0.25, y = -.25)) %>%
        plotly::style(hoverinfo = 'none')
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
                    # shiny::checkboxInput("disc", "Discontinuous y-axis", FALSE),
                    shiny::numericInput("thresh", "Hide between threshold", 0, step = 0.01),
                    shiny::sliderInput("Xrange", "X-range:", min = -10, max = 10, value = c(-5, 5)),
                    shiny::sliderInput("numrow",label =  "Number of rows:", min = 1, max = 10, value = 1),
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
                    plotly::plotlyOutput("distPlotly"))
    )
  )


  shiny::shinyApp(ui = ui, server = server)
}
