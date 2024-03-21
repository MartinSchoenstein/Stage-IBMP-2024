PATH = "/Users/jroignant/Desktop/Stage_Martin/complete_app/RNAvis/data/data_reduice/FLEPseq_runs"
listdir = list.dirs(path = PATH, full.names = FALSE ,recursive = FALSE)
barcodeFileSimilar = read_tsv("/Users/jroignant/Desktop/Stage_Martin/complete_app/RNAvis/data/data_reduice/FLEPseq_runs/02_RUN02_20221021/4_Tail/barcode08_reduice.read_info.result.merged.parts.csv")

ui <- dashboardPage(
  dashboardHeader(disable = TRUE
  ),
  dashboardSidebar(
    sidebarMenu(id = "sidebar",
      menuItem("  Single Gene", tabName = "singleGene", icon = icon("play")) ,
      menuItem("  Multiple Genes", tabName = "multipleGene", icon = icon("forward")) ,
      menuItem("  Similarities", tabName = "similarGene", icon = icon("spinner")),
      textInput("yMax", "Y Range", value = 0.03)
  )),
  dashboardBody(
    tabItems(
      tabItem(tabName = "singleGene",
        selectInput("choixRun", "Selectionnez un run :", listdir),
        selectInput("choixGene", "Selectionnez un gène :", c()),
        plotOutput("plotPolya"),
        checkboxInput("meanPlot", "Mean Plot", value = FALSE),
        downloadButton("downloadPlotPolya", "Download this plot")
      ),
      
      tabItem(tabName = "multipleGene",
        selectInput("choixRunMultiple", "Selectionnez un run :", listdir),
        selectInput("choixGeneMultiple", "Selectionnez des gènes :", c(), multiple = TRUE),
        plotOutput("plotPolyaMultiple"),
        downloadButton("downloadPlotPolyaMultiple", "Download this plot")
      ),
      
      tabItem(tabName = "similarGene",
        p("Choisissez votre profil de référence :"),
        selectInput("choixRunSimilar", "Selectionnez un run :", listdir),
        selectInput("choixBarcodeSimilar", "Selectionnez un barcode :", c(), width = 360),
        selectInput("choixGeneSimilar", "Selectionnez un gène : ", c()),
        plotOutput("plotReferenceSimilar"),
        p(""),
        actionButton("buttonReference", label = "Valider la référence", icon = icon("check")),
        actionButton("buttonSimilar", label = "Lancer la recherche de gènes similaires", icon = icon("magnifying-glass-chart")),
        textOutput("meanReference"),
        textOutput("csvSimilar")
      )
  ))
)

server <- function(input, output) {
  
  #----------------------------------------Single Gene
  
  pathBarcodeReact = reactive({
   pathBarcode =  paste0(PATH, "/", input$choixRun, "/4_Tail")
  })
  
  listCsvReact = reactive({
    listCsv = list.files(pathBarcodeReact(), full.names = FALSE, pattern = "csv$")
  })
  
  bindBarcodeFileReact = reactive({
    bindBarcodeFile = rbindlist(lapply(paste0(pathBarcodeReact(), "/", listCsvReact()), fread), idcol = "barcode")
  })
  
  observeEvent(input$choixRun, {
    listGene = unique(bindBarcodeFileReact()$mRNA)
    updateSelectInput(inputId = "choixGene", choices = listGene)
  })
  
  subdfBindBarcodeFileReact = reactive ({
    subdfBindBarcodeFile = subset.data.frame(bindBarcodeFileReact(), mRNA == input$choixGene)
  })
  
  observe({
    if (input$meanPlot){
      output$plotPolya = renderPlot({
        print(plotPolyaMeanReact())
          })
    }
    else{
      output$plotPolya = renderPlot({
        print(plotPolyaReact())
      })}
    
    })
  
  plotPolyaReact = reactive({
    plotPolya = ggplot(subdfBindBarcodeFileReact(), aes(x = polya_length, group = barcode, color = as.character(barcode))) + geom_density() + xlim(0, 250) + ylim(0, as.numeric(input$yMax)) +
      scale_color_hue(name = "", labels = paste(listCsvReact())) +
      xlab("PolyA tail length") + ylab("Density") +
      theme(legend.title = element_text(size = 14), legend.text = element_text(size = 10), legend.position = "bottom")
  })
  
  plotPolyaMeanReact = reactive({
    plotPolyaMean = ggplot(subdfBindBarcodeFileReact(), aes(x = polya_length)) + geom_density() + xlim(0, 250) + ylim(0, as.numeric(input$yMax)) +
      xlab("PolyA tail length") + ylab("Density") +
      labs(caption = "") +
      theme(plot.caption = element_text(size = 27 + 4.8*length(listCsvReact())))
  })
  
  output$downloadPlotPolya = downloadHandler(
    filename = function(){
      paste0("plot-",input$choixRun, "-", input$choixGene, ".png")
      },
    content = function(file){
      if (input$meanPlot){
        ggsave(file, plotPolyaMeanReact(), width = 16)}
      else {
        ggsave(file, plotPolyaReact(), width = 16)
      }
  })

#----------------------------------------Multiple Gene  
  
  pathBarcodeMultipleReact = reactive({
    pathBarcodeMultiple =  paste0(PATH, "/", input$choixRunMultiple, "/4_Tail")
  })
  
  listCsvMultipleReact = reactive({
    listCsvMultiple = list.files(pathBarcodeMultipleReact(), full.names = FALSE, pattern = "csv$")
  })
  
  bindBarcodeFileMultipleReact = reactive({
    bindBarcodeFileMultiple = rbindlist(lapply(paste0(pathBarcodeMultipleReact(), "/", listCsvMultipleReact()), fread), idcol = "barcode")
  })
  
  observeEvent(input$choixRunMultiple, {
    listGeneMultiple = unique(bindBarcodeFileMultipleReact()$mRNA)
    updateSelectInput(inputId = "choixGeneMultiple", choices = listGeneMultiple)
  })
  
  subdfBindBarcodeFileMultipleReact = reactive ({
    subdfBindBarcodeFileMultiple = subset.data.frame(bindBarcodeFileMultipleReact(), mRNA %in% input$choixGeneMultiple)
  })
  
  plotPolyaMultipleReact = reactive({
    plotPolya = ggplot(subdfBindBarcodeFileMultipleReact(), aes(x = polya_length, color = mRNA)) + geom_density() + xlim(0, 250) + ylim(0, as.numeric(input$yMax)) +
      scale_color_hue(name = "", labels = input$choixGeneMultiple) +
      xlab("PolyA tail length") + ylab("Density") +
      theme(legend.title = element_text(size = 14), legend.text = element_text(size = 10), legend.position = "bottom") +
      facet_wrap(vars(barcode))
  })
  
  observe({
      output$plotPolyaMultiple = renderPlot({
        print(plotPolyaMultipleReact())
      })
  })
  
  output$downloadPlotPolyaMultiple = downloadHandler(
    filename = function(){
      paste0("plotMultiple-",input$choixRunMultiple, "-", paste(input$choixGeneMultiple, collapse = "-"), ".png")
    },
    content = function(file){
        ggsave(file, plotPolyaMultipleReact(), width = 16)
    })
  
#-------------------------------------Similar
  
  pathBarcodeSimilarReact = reactive({
    pathBarcodeSimilar =  paste0(PATH, "/", input$choixRunSimilar, "/4_Tail")
  })
  
  listCsvSimilarReact = reactive({
    listCsvSimilar = list.files(pathBarcodeSimilarReact(), full.names = FALSE, pattern = "csv$")
  })
  
  observeEvent(input$choixRunSimilar, {
    updateSelectInput(inputId = "choixBarcodeSimilar", choices = listCsvSimilarReact())
  })
  
  barcodeFileSimilarReact = reactive({
    barcodeFileSimilar = read_tsv(paste0(pathBarcodeSimilarReact(), "/", input$choixBarcodeSimilar))
  })
  
  subdfBarcodeFileSimilarReact = reactive ({
    subdfBarcodeFileSimilar = subset.data.frame(barcodeFileSimilarReact(), mRNA == input$choixGeneSimilar)
  })
  
  observeEvent(input$choixBarcodeSimilar, {
    listGeneSimilar = unique(barcodeFileSimilar$mRNA)
    updateSelectInput(inputId = "choixGeneSimilar", choices = listGeneSimilar)
  })
  
  plotPolyaSimilarReact = reactive({
    plotPolyaSimilar = ggplot(subdfBarcodeFileSimilarReact(), aes(x = polya_length)) + geom_density() + xlim(0, 250) + ylim(0, as.numeric(input$yMax)) +
      xlab("PolyA tail length") + ylab("Density")
  })
  
  output$plotReferenceSimilar = renderPlot({
    print(plotPolyaSimilarReact())
  })
  
  observeEvent(input$buttonReference, {
    listCsvReference = listCsvSimilarReact()
    fileReference = subdfBarcodeFileSimilarReact()
    plotReference = plotPolyaSimilarReact()
    meanReference = max(plotPolyaSimilarReact())   #marche pas 
    output$meanReference = renderText({
      print(meanReference)
        })
  })
  
  observeEvent(input$buttonSimilar, {  #ajouter une manière de contrôler que la reference a bien été mise (si on reste sur l'option de 2 boutons) sinon si on clique d'abord à droite ça fait crash l'appli
    for(a in listCsvReference){           #marche si on met listCsvSimilarReact() direct mais ducoup si je bouge les input après avoir validé ça passe au dessus de la validation
      output$csvSimilar = renderText({
        print(a)
      })
    }
  })
  
}
shinyApp(ui, server)


