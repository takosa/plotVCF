library(shiny)
library(VariantAnnotation)
library(karyoploteR)
library(ggplot2)
library(magrittr)
#library(plotly)

ui <- fluidPage(

    titlePanel("plot VCF"),

    sidebarLayout(
        sidebarPanel(width = 2,
            fileInput("inFile", "Input VCF", accept = c("vcf", "vcf.gz", "bcf", "bcf.gz"))
        ),
        mainPanel(width = 10,
            tabsetPanel(
                tabPanel("genotype", 
                    plotOutput("geno", height = "800px")
                ),
                tabPanel("position", 
                    plotOutput("pos", height = "800px")
                ),
                tabPanel("MAF", 
                    plotOutput("maf", height = "800px")
                ),
                tabPanel("yield", 
                    tableOutput("summary"),
                    plotOutput("yield", height = "800px")
                ),
                tabPanel("hint", 
                    includeMarkdown("hint.md")
                )
            )
        )
    )
)

server <- function(input, output) {
    
    readFile <- reactive({
        req(input$inFile)
        vcf <- readVcf(input$inFile$datapath)
        vcf
    })
    
    output$geno <- renderPlot({
        vcf <- readFile()
        gt <- geno(vcf)$GT
        gt <- reshape2::melt(gt, varnames = c("marker", "sample"), value.name="genotype")
        gp <- ggplot(gt, aes(x = sample, y = marker)) + 
            geom_raster(aes(fill = genotype)) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
        gp
    })
    
    output$pos <- renderPlot({
        vcf <- readFile()
        contig <- meta(header(vcf))$contig
        genome <- toGRanges(data.frame(
            chr = rownames(contig),
            start = 1L,
            end = contig$length
        ))
        params <- getDefaultPlotParams(plot.type = 1)
        params$leftmargin <- 0.2
        plotKaryotype(genome = genome, plot.type = 1, plot.params = params) %>% 
            kpAddBaseNumbers() %>% 
            kpPlotMarkers(data = rowRanges(vcf), labels = rowRanges(vcf)$paramRangeID)
    
    })
    output$maf <- renderPlot({
        vcf <- readFile()
        info <- info(vcf)
        validate(
            need(!is.null(input$MAF), "There is no INFO/MAF tag. To run `bcftools plugin fill-tags -o out.vcf in.vcf` may be helpful.")
        )
        d <- data.frame(marker = rownames(vcf), MAF = unlist(info$MAF))
        ggplot(d, aes(x = marker, y = MAF)) +
            geom_point() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    })
    
    output$summary <- renderTable({
        vcf <- readFile()
        gt <- geno(vcf)$GT
        all_datapoint <- as.integer(prod(dim(vcf)))
        called_datapoint <- sum(gt != "./.")
        data.frame("All data points" = all_datapoint,
                   "Called data points" = called_datapoint,
                   "Yield(%)" = called_datapoint/all_datapoint, 
                   check.names = FALSE)
    })
    output$yield <- renderPlot({
        vcf <- readFile()
        gt <- geno(vcf)$GT
        ratio_marker <- apply(gt, 1, function(x)  sum(x != "./.") / length(x) )
        ratio_sample <- apply(gt, 2, function(x)  sum(x != "./.") / length(x) )
        type <- rep(c("marker-wise", "sample-wise"), c(length(ratio_marker), length(ratio_sample)))
        name <- c(rownames(vcf), colnames(vcf))
        fraction_no_missing <- c(ratio_marker, ratio_sample)
        d <- data.frame(type = type, name = name, yield = fraction_no_missing)
        ggplot(d, aes(name, yield)) +
            geom_point() +
            facet_wrap(vars(type), ncol = 1, scales = "free_x") +
            geom_hline(aes(yintercept = mean(yield))) +
            ylim(0, 1) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    })

}

shinyApp(ui = ui, server = server)
