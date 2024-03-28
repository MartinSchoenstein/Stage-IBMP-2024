PATH = "/Users/jroignant/Desktop/Stage_Martin/complete_app/RNAvis/data/data_reduice/FLEPseq_runs"
listdir = list.dirs(path = PATH, full.names = FALSE ,recursive = FALSE)
PATHo = "/Users/jroignant/Desktop/Stage_Martin/complete_app/RNAvis/data/data/FLEPseq_runs"

ui <- dashboardPage(skin = "purple",
  
  dashboardHeader(title = "Project"
    
  ),
#========= Sidebar ===========
  dashboardSidebar(
    sidebarMenu(id = "sidebar",
      menuItem("", tabName = "home", icon = icon("home", style = "color: #9d52ff")),
      menuItem("  .Single Gene Overview", tabName = "single", icon = icon("play", style = "color: #9d52ff")) ,
      menuItem("  .Multiple Genes Overview", tabName = "multiple", icon = icon("forward", style = "color: #9d52ff")) ,
      menuItem("  .CCR4 Mutation", tabName = "mutation", icon = icon("seedling", style = "color: #9d52ff"), badgeLabel = " ", badgeColor = "red"),
      conditionalPanel(condition = "input.sidebar == 'single'",
          selectInput("runSingle", "Choisissez un run :", listdir),
          selectInput("geneSingle", "Choisissez un gène :", c()),
          numericInput("seuilPercU", "Minimal percentage of U :", value = 60),
          actionButton("buttonValidationSingle", label = "Valider cette selection", icon = icon("check", style = "color: #1dcc1d"))
      ),
      conditionalPanel(condition = "input.sidebar == 'multiple'",
          selectInput("runMultiple", "Choisssez un run", listdir),
          numericInput("seuilPercUMultiple", "Minimal percentage of U :", value = 60),
          actionButton("buttonValidationMultiple", label = "Valider cette selection", icon = icon("check", style = "color: #1dcc1d"))
      )
  )),
#====== Body ======= 
  dashboardBody(
    tabItems(
      tabItem(tabName = "single",
        h1("Single Gene Overview"),
        tableOutput("tableStatReads"),      
        #h3("Uridylation Statistics :"),
        #p(""),
        box(title = "Percentage of Uridylation :", status = "primary", solidHeader = TRUE, width = 12 , plotOutput("uridylationPercentageSingle"), downloadButton("downloadUrdidylPercPlot", "Download this plot")),
        #checkboxInput("buttonChangePlotUrdidyl", "Change plot style", value = FALSE)
        box(title = "Number of U in the added tail :",  status = "primary", solidHeader = TRUE, width = 12, plotOutput("uridylationNumberSingle"), downloadButton("downloadUrdidylNumberPlot", "Download this plot")),
        #h3("Poly(A) tail :"),
        #p(""),
        box(title = "Poly(A) tail length distribution :", status = "primary", solidHeader = TRUE, width = 12, plotOutput("polyaDistributionSingle"), radioButtons("plotType", "Change plot", c("Density plot", "Cumulativ plot", "Histogram plot"), selected = "Density plot"), downloadButton("downloadPolyaDistribPlot", "Download this plot")),
        #p(""),
        box(title = "Poly(A) tail length distribution in uridylate mRNA", status = "primary", solidHeader = TRUE, width = 12, plotOutput("polyaDistributionUrdiylSingle"), radioButtons("plotTypeUridyl", "Change plot", c("Density plot", "Cumulativ plot", "Histogram plot"), selected = "Density plot"), downloadButton("downloadPolyaDistributionUridylPlot", "Download this plot")),
        p(".")
      ),
      
      tabItem(tabName = "multiple",
        h1("Multiple Gene Overview"),
        tableOutput("tableStatReadsMultiple"),
        box(title = "Percentage of Uridylation :", status = "primary", solidHeader = TRUE, width = 12, plotOutput("uridylationPercentageMultiple"), downloadButton("downloadUridylPercMultiplePlot", "Download this plot")),
        box(title = "Number of U in the added tail :", status = "primary", solidHeader = TRUE, width = 12, plotOutput("uridylationNumberMultiple"), downloadButton("downloadUridylNumberMultiplePlot", "Download this plot")),
        box(title = "Poly(A) tail length distribution :", status = "primary", solidHeader = TRUE, width = 12, plotOutput("polyaDistributionMultiple"), radioButtons("plotTypeMultiple", "Change plot", c("Density plot", "Cumulativ plot", "Histogram plot"), selected = "Density plot"), downloadButton("downloadPolyaDistribMultiplePlot", "Download this plot")),
        box(title = "Poly(A) tail length distribution in uridylate mRNA", status = "primary", solidHeader = TRUE, width = 12, plotOutput("polyaDistributionUridylMultiple"), radioButtons("plotTypeUridylMultiple", "Change plot", c("Density plot", "Cumulativ plot", "Histogram plot"), selected = "Density plot"), downloadButton("downloadPolyaDistributionUridylMultiplePlot", "Download this plot")),
        p(".")
      )
    )
  )
)

server <- function(input, output) { 

#========================= Mise en place ===========================
    
  listCsvReact = reactive({                #liste réatcive des différents fichier csv pour le run choisi
    listCsv = list.files(paste0(PATH, "/", input$runSingle, "/4_Tail"), full.names = FALSE, pattern = "csv$")
  })
  
  listBarcodeReact = eventReactive(input$buttonValidationSingle, {            #liste réactive des différents barcodes présents
    listBarcode = c()
    for (b in listCsvReact()){
      listBarcode = append(listBarcode, substring(b, 1, 9))
    }
    return (listBarcode)
  })
  
  bindCsvFilesReact = reactive({           #fusion en un tableau réactif des fichiers csv et ajout d'une ligne barcode (valeure 1 à x en fonction de l'ordre alphabétique)
    bindCsvFiles = rbindlist(lapply(paste0(paste0(PATH, "/", input$runSingle, "/4_Tail"), "/", listCsvReact()), fread), idcol = "barcode")
  })

  dfUridylReact = eventReactive(input$buttonValidationSingle, {      #tableau modifié en ajoutant une colonne pour différencier les lignes uridylés des autres (selon le pourcentage d'input) et on ne garde que le gène séléctionné, c'est ce tableau qu'on utilisera après, ne se met à jour que lors de la validation des inputs
    dfUridyl = bindCsvFilesReact()%>%mutate(uridyl = if_else(add_tail_pct_T>input$seuilPercU, "U", "/"))%>%filter(mRNA == input$geneSingle)
  })
  
  observeEvent(input$runSingle, {          #mise à jour du choix de gène en fonction de ceux présents dans le fichier fusionné
    updateSelectInput(inputId = "geneSingle", choices = unique(bindCsvFilesReact()$mRNA))
  })
  
  # - - - - - - - - - 
  
  listCsvMultipleReact = reactive({                #liste réatcive des différents fichier csv pour le run choisi
    listCsvmultiple = list.files(paste0(PATH, "/", input$runMultiple, "/4_Tail"), full.names = FALSE, pattern = "csv$")
  })
  
  listBarcodeMultipleReact = eventReactive(input$buttonValidationMultiple, {            #liste réactive des différents barcodes présents qui se met à jour lors de la validation
    listBarcodeMultiple = c()
    for (b in listCsvMultipleReact()){
      listBarcodeMultiple = append(listBarcodeMultiple, substring(b, 1, 9))
    }
    return (listBarcodeMultiple)
  })
  
  bindCsvFilesMultipleReact = reactive({           #fusion en un tableau réactif des fichiers csv et ajout d'une ligne barcode (valeure 1 à x en fonction de l'ordre alphabétique)
    bindCsvFilesMultiple = rbindlist(lapply(paste0(paste0(PATH, "/", input$runMultiple, "/4_Tail"), "/", listCsvMultipleReact()), fread), idcol = "barcode")
  })
  
  dfUridylMultipleReact = eventReactive(input$buttonValidationMultiple, {      #tableau modifié en ajoutant une colonne pour différencier les lignes uridylés des autres (selon le pourcentage d'input) et mais on garde tout le gène cette fois ci, c'est ce tableau qu'on utilisera après pour le multiple, ne se met à jour que lors de la validation des inputs
    dfUridylMultiple = bindCsvFilesMultipleReact()%>%mutate(uridyl = if_else(add_tail_pct_T>input$seuilPercUMultiple, "U", "/"))%>%mutate(x = ".")
  })
  
  
  #========================================================= Single =================================================================
  
  
  #======================== Bouton Valider ============================
  
  observeEvent(input$buttonValidationSingle, {     #tout ce qui se passe qu'une fois que la séléctionner du run et gène est validé
  
    # ---------->Tableau récap pour le run (barcode et génotype) et le gène (nombre de reads) choisis
    
    listBarcodeCorrespondance = list.files(paste0(PATHo, "/", input$runSingle),full.names = TRUE, pattern = "^barcode_correspondance")
    barcodeCorrespondance = read_tsv(listBarcodeCorrespondance, col_names = FALSE)%>%arrange(X1)     #on trie pour avoir les barcodes dans le même sens alphabétique dans notre tableau
    ReadsNumber = c()      #nombre de reads totals de ce gène
    ReadsUridylNumber = c()  #nombre de reads totals de ce gène qui sont considérés uridylés
    for (b in 1:nrow(barcodeCorrespondance)){
      ReadsNumber = append(ReadsNumber, nrow(dfUridylReact()%>%filter(barcode == b)))
      ReadsUridylNumber = append(ReadsUridylNumber, nrow(dfUridylReact()%>%filter(barcode == b)%>%filter(uridyl == "U")))
    }
    tableStatReads = barcodeCorrespondance%>%mutate(ReadsNumber, ReadsUridylNumber)
    if (nrow(tableStatReads) == 3){         #si le tableau n'a que trois colonne ça veut dire que dans le fichier de base barcode_correspondance on a le barcore et le genotype dans la même colonne donc on sépare
      tableStatReads = tableStatReads%>%separate(X1, c("Barcode", "Genotype"))
      colnames(tableStatReads)[3] = "Number of reads"
      colnames(tableStatReads)[4] = "Number of reads with uridylation"
    }
    else{  #sinon on les renomme juste
      colnames(tableStatReads)[1] = "Barcode"
      colnames(tableStatReads)[2] = "Genotype"
      colnames(tableStatReads)[3] = "Number of reads"
      colnames(tableStatReads)[4] = "Number of reads with uridylation"
    }
    
    # ---------->Tous les outputs
    
    output$tableStatReads = renderTable(tableStatReads)       #table des stats
    
    output$uridylationPercentageSingle = renderPlot({         #plot du pourcentage d'uridylation par génotype pour un gène choisis
      print(plotUridylationPercentageSingleReact())
      })
    
    output$uridylationNumberSingle = renderPlot({             #plot du nombre de U ajouté par génotype pour un gène choisis
      print(plotUridylationNumberSingleReact())
      })
    
    observe({                                                 #le plot dépend de l'input présent sous le graphique, on affiche le densité par défault
      if (input$plotType == "Density plot"){
        output$polyaDistributionSingle = renderPlot({
          print(plotPolyaDistributionSingleReact())
      })
      }
      if (input$plotType == "Cumulativ plot"){
        output$polyaDistributionSingle = renderPlot({
          print(plotPolyaDistributionSingle2React())
        })
      }
      if (input$plotType == "Histogram plot"){
        output$polyaDistributionSingle = renderPlot({
          print(plotPolyaDistributionSingle3React())
        })
      }
      })
    
    observe({                                               #pareil
      if (input$plotTypeUridyl == "Density plot"){
        output$polyaDistributionUrdiylSingle = renderPlot({
          print(plotPolyaDistributionUridylSingleReact())
          })
        }
      if (input$plotTypeUridyl == "Cumulativ plot"){
        output$polyaDistributionUrdiylSingle = renderPlot({
          print(plotPolyaDistributionUridylSingle2React())
          })
        }
      if (input$plotTypeUridyl == "Histogram plot"){
        output$polyaDistributionUrdiylSingle = renderPlot({
          print(plotPolyaDistributionUridylSingle3React())
          })
        }
      })
  })

  #===================== Plot Pourcentage Uridylation ===========================
  
  dfPercUridylReact = reactive({      #création d'un tableau qui servira à construire le graph, on calcule notamment le pourcentage de "lignes uridylées", la colonne position permettra de placer le texte selon la valeur dans le graphique
    mRNA = c()
    Uridylation = c()
    status = c()
    barcode = c()
    position = c()
    for (br in unique(dfUridylReact()$barcode)){
      mRNA = append(mRNA, input$geneSingle)
      mRNA = append(mRNA, input$geneSingle)
      Uridylation = append(Uridylation, nrow(dfUridylReact()%>%filter(uridyl == "/")%>%filter(barcode == br)) / nrow(dfUridylReact()%>%filter(barcode == br)) * 100)
      Uridylation = append(Uridylation, nrow(dfUridylReact()%>%filter(uridyl == "U")%>%filter(barcode == br)) / nrow(dfUridylReact()%>%filter(barcode == br)) * 100)
      status = append(status, "no Uridylation")
      status = append(status, "Uridylation")
      barcode = append(barcode, br)
      barcode = append(barcode, br)
      position = append(position, nrow(dfUridylReact()%>%filter(uridyl == "/")%>%filter(barcode == br)) / nrow(dfUridylReact()%>%filter(barcode == br)) * 100 + nrow(dfUridylReact()%>%filter(uridyl == "U")%>%filter(barcode == br)) / nrow(dfUridylReact()%>%filter(barcode == br)) * 100)
      position = append(position, nrow(dfUridylReact()%>%filter(uridyl == "U")%>%filter(barcode == br)) / nrow(dfUridylReact()%>%filter(barcode == br)) * 100)
      }
    dfPercUridyl = data.frame(mRNA, Uridylation, status, barcode, position)
    return(dfPercUridyl)
  })
  
  plotUridylationPercentageSingleReact = reactive({           #Plot pourcentage Uridylation
    plotUridylationPercentageSingle = ggplot(data = dfPercUridylReact()) +
      geom_bar(aes(x = mRNA, y = Uridylation, color = as.character(barcode), fill = status), stat = "identity") + 
      facet_wrap(vars(barcode)) + 
      scale_color_hue(name = "", label = paste(listBarcodeReact())) +
      xlab("Gene (split in different barcodes)") + ylab("Percentage of presence (%)") +
      scale_fill_manual(name = "Additional Tail", values = c("#f1e3ca", "#87E990")) +
      ylim(0, 50) +
      theme(legend.position="bottom", legend.box = "vertical") +
      geom_text(x = 1, y = dfPercUridylReact()$position, label = paste(round(dfPercUridylReact()$Uridylation, 2), "%"), size = (5 * (dfPercUridylReact()$Uridylation/(dfPercUridylReact()$Uridylation - 0.005))), vjust = 1.05) #le calcul pour la taille permet d'avoo-ir une taille nulle lorsqu'on a 0% mais influe presque pas sur la taille sinon
  })

  output$downloadUrdidylPercPlot = downloadHandler(        #Télechargement du plot en image
    filename = function(){
      paste0("plotPercentageUridylation-",input$runSingle, "-", input$geneSingle, ".png")
    },
    content = function(file){
      ggsave(file, plotUridylationPercentageSingleReact(), width = 16)
    })
  
  #===================== Plot Nombre U ============================
  
  plotUridylationNumberSingleReact = reactive({               #Plot nombre de U ajouté 
    plotUridylationNumberSingle = ggplot(data = dfUridylReact(), aes(x = mRNA, y = add_tail_T, color = as.character(barcode))) +
      geom_violin(trim = FALSE) +
      facet_wrap(vars(barcode)) +
      geom_boxplot(width = 0.15) +
      ylim(0,10) +
      scale_color_hue(name = "", label = paste(listBarcodeReact())) +
      theme(legend.position="bottom", legend.box = "vertical") +
      xlab("Gene (split in different barcodes)") + ylab("Number of U")
  })
    
  output$downloadUrdidylNumberPlot = downloadHandler(        #Télechargement du plot en image
    filename = function(){
      paste0("plotNumberUridyl-",input$runSingle, "-", input$geneSingle, ".png")
    },
    content = function(file){
      ggsave(file, plotUridylationNumberSingleReact(), width = 16)
    })
  
  # =========================== Plot Distribution taille queue polyA ==============================
  

    plotPolyaDistributionSingleReact = reactive({                 #Plot taille de la queue polyA (densité)
      plotPolyaDistributionSingle = ggplot(dfUridylReact(), aes(x = polya_length, group = barcode, color = as.character(barcode))) + geom_density() + xlim(10, 200) +
        scale_color_hue(name = "", labels = paste(listBarcodeReact())) +
        xlab("PolyA tail length") + ylab("Density of distribution") +
        theme(legend.title = element_text(size = 14), legend.text = element_text(size = 10), legend.position = "bottom")
    })
  
    plotPolyaDistributionSingle2React = reactive({                #Plot taille de la queue polyA (cumulatif)
      plotPolyaDistributionSingle = ggplot(dfUridylReact(), aes(polya_length, color = as.character(barcode))) + stat_ecdf(geom = "line") +
        scale_color_hue(name = "", labels = paste(listBarcodeReact())) +
        theme(legend.title = element_text(size = 14), legend.text = element_text(size = 10), legend.position = "bottom") +
        xlab("PolyA tail length") + ylab("Cumulative proportion") 
    })
    
    plotPolyaDistributionSingle3React = reactive({                #Plot taille de la queue polyA (boxplot)
      plotPolyaDistributionSingle = ggplot(data = dfUridylReact(), aes(x = mRNA, y = polya_length, color = as.character(barcode))) +
        geom_violin(trim = FALSE) +
        facet_wrap(vars(barcode)) +
        geom_boxplot(width = 0.15) +
        ylim(0,200) +
        scale_color_hue(name = "", labels = paste(listBarcodeReact())) +
        theme(legend.position="bottom") +
        xlab("Gene (split in different barcodes)") + ylab("Poly(A) tail length")
    })
    
    output$downloadPolyaDistribPlot = downloadHandler(           #Télechargement du plot en image, le plot téléchargé dépend de celui affiché donc de la valeur de l'input
      filename = function(){
        paste0("plotDistributionPolya-",input$runSingle, "-", input$geneSingle, ".png")
      },
      content = function(file){
        if (input$plotType == "Density plot"){
          ggsave(file, plotPolyaDistributionSingleReact(), width = 12)}
        else if (input$plotType == "Cumulativ plot"){
          ggsave(file, plotPolyaDistributionSingle2React(), width = 12)}
        else {
          ggsave(file, plotPolyaDistributionSingle3React(), width = 12)}
      })

# ===================== Plot Distribution des tailles de queues polyA des mRNA uridylés ==============    
    
    plotPolyaDistributionUridylSingleReact = reactive({                #Plot taille de la queue polyA (densité)
      polyaDistributionUridylSingle = ggplot(dfUridylReact()%>%filter(uridyl == "U"), aes(x = polya_length, group = barcode, color = as.character(barcode))) + geom_density() + xlim(0, 200) +
        scale_color_hue(name = "", labels = paste(listBarcodeReact())) +
        xlab("PolyA tail length") + ylab("Density of distribution") +
        theme(legend.title = element_text(size = 14), legend.text = element_text(size = 10), legend.position = "bottom")
    })
 
    plotPolyaDistributionUridylSingle2React = reactive({               #Plot taille de la queue polyA (cumulatif)
      plotPolyaDistributionUridylSingle = ggplot(dfUridylReact()%>%filter(uridyl == "U"), aes(polya_length, color = as.character(barcode))) + stat_ecdf(geom = "line") +
        scale_color_hue(name = "", labels = paste(listBarcodeReact())) +
        theme(legend.title = element_text(size = 14), legend.text = element_text(size = 10), legend.position = "bottom") +
        xlab("PolyA tail length") + ylab("Cumulative proportion") 
    })
       
    plotPolyaDistributionUridylSingle3React = reactive({               #Plot taille de la queue polyA (boxplot)
      plotPolyaDistributionUridylSingle = ggplot(data = dfUridylReact()%>%filter(uridyl == "U"), aes(x = mRNA, y = polya_length, color = as.character(barcode))) +
        geom_violin(trim = FALSE) +
        facet_wrap(vars(barcode)) +
        geom_boxplot(width = 0.15) +
        ylim(0,250) +
        scale_color_hue(name = "", labels = paste(listBarcodeReact())) +
        theme(legend.position="bottom") +
        xlab("Gene (split in different barcodes)") + ylab("Poly(A) tail length")
    })
    
    
#========================================================= Multiple =================================================================    

    #======================== Bouton Valider ============================
    
    observeEvent(input$buttonValidationMultiple, {     #tout ce qui se passe qu'une fois que la séléctionner du run est validé
      
      # ---------->Tableau récap pour le run (barcode et génotype, nombre de reads)
      
      listBarcodeCorrespondanceMultiple = list.files(paste0(PATHo, "/", input$runMultiple),full.names = TRUE, pattern = "^barcode_correspondance")
      barcodeCorrespondanceMultiple = read_tsv(listBarcodeCorrespondanceMultiple, col_names = FALSE)%>%arrange(X1)     #on trie pour avoir les barcodes dans le même sens alphabétique dans notre tableau
      ReadsNumberMultiple = c()      #nombre de reads totals 
      ReadsUridylNumberMultiple = c()  #nombre de reads totals qui sont considérés uridylés
      for (a in 1:nrow(barcodeCorrespondanceMultiple)){
        ReadsNumberMultiple = append(ReadsNumberMultiple, nrow(dfUridylMultipleReact()%>%filter(barcode == a)))
        ReadsUridylNumberMultiple = append(ReadsUridylNumberMultiple, nrow(dfUridylMultipleReact()%>%filter(barcode == a)%>%filter(uridyl == "U")))
      }
      tableStatReadsMultiple = barcodeCorrespondanceMultiple%>%mutate(ReadsNumberMultiple, ReadsUridylNumberMultiple)
      if (nrow(tableStatReadsMultiple) == 3){         #si le tableau n'a que trois colonne ça veut dire que dans le fichier de base barcode_correspondance on a le barcore et le genotype dans la même colonne donc on sépare
        tableStatReadsMultiple = tableStatReadsMultiple%>%separate(X1, c("Barcode", "Genotype"))
        colnames(tableStatReadsMultiple)[3] = "Number of reads"
        colnames(tableStatReadsMultiple)[4] = "Number of reads with uridylation"
      }
      else{  #sinon on les renomme juste
        colnames(tableStatReadsMultiple)[1] = "Barcode"
        colnames(tableStatReadsMultiple)[2] = "Genotype"
        colnames(tableStatReadsMultiple)[3] = "Number of reads"
        colnames(tableStatReadsMultiple)[4] = "Number of reads with uridylation"
      }
      
      # ---------->Tous les outputs
      
      output$tableStatReadsMultiple = renderTable(tableStatReadsMultiple)       #table des stats
      
      output$uridylationPercentageMultiple = renderPlot({         #plot du pourcentage d'uridylation par génotype 
        print(plotUridylationPercentageMultipleReact())
      })
      
      output$uridylationNumberMultiple = renderPlot({             #plot du nombre de U ajouté par génotype 
        print(plotUridylationNumberMultipleReact())
      })
      
      observe({                                                 #le plot dépend de l'input présent sous le graphique, on affiche le densité par défault
        if (input$plotTypeMultiple == "Density plot"){
          output$polyaDistributionMultiple = renderPlot({
            print(plotPolyaDistributionMultipleReact())
          })
        }
        if (input$plotTypeMultiple == "Cumulativ plot"){
          output$polyaDistributionMultiple = renderPlot({
            print(plotPolyaDistributionMultiple2React())
          })
        }
        if (input$plotTypeMultiple == "Histogram plot"){
          output$polyaDistributionMultiple = renderPlot({
            print(plotPolyaDistributionMultiple3React())
          })
        }
      })
      
      observe({                                               #pareil
        if (input$plotTypeUridylMultiple == "Density plot"){
          output$polyaDistributionUridylMultiple = renderPlot({
            print(plotPolyaDistributionUridylMultipleReact())
          })
        }
        if (input$plotTypeUridylMultiple == "Cumulativ plot"){
          output$polyaDistributionUridylMultiple = renderPlot({
            print(plotPolyaDistributionUridylMultiple2React())
          })
        }
        if (input$plotTypeUridylMultiple == "Histogram plot"){
          output$polyaDistributionUridylMultiple = renderPlot({
            print(plotPolyaDistributionUridylMultiple3React())
          })
        }
      })
    })
    
    #===================== Plot Pourcentage Uridylation ===========================
    
    dfPercUridylMultipleReact = reactive({      #création d'un tableau qui servira à construire le graph, on calcule notamment le pourcentage de "lignes uridylées", la colonne position permettra de placer le texte selon la valeur dans le graphique
      x = c()
      Uridylation = c()
      status = c()
      barcode = c()
      position = c()
      for (br in unique(dfUridylMultipleReact()$barcode)){
        x = append(x, ".")
        x = append(x, ".")
        Uridylation = append(Uridylation, nrow(dfUridylMultipleReact()%>%filter(uridyl == "/")%>%filter(barcode == br)) / nrow(dfUridylMultipleReact()%>%filter(barcode == br)) * 100)
        Uridylation = append(Uridylation, nrow(dfUridylMultipleReact()%>%filter(uridyl == "U")%>%filter(barcode == br)) / nrow(dfUridylMultipleReact()%>%filter(barcode == br)) * 100)
        status = append(status, "no Uridylation")
        status = append(status, "Uridylation")
        barcode = append(barcode, br)
        barcode = append(barcode, br)
        position = append(position, nrow(dfUridylMultipleReact()%>%filter(uridyl == "/")%>%filter(barcode == br)) / nrow(dfUridylMultipleReact()%>%filter(barcode == br)) * 100 + nrow(dfUridylMultipleReact()%>%filter(uridyl == "U")%>%filter(barcode == br)) / nrow(dfUridylMultipleReact()%>%filter(barcode == br)) * 100)
        position = append(position, nrow(dfUridylMultipleReact()%>%filter(uridyl == "U")%>%filter(barcode == br)) / nrow(dfUridylMultipleReact()%>%filter(barcode == br)) * 100)
      }
      dfPercUridylMultiple = data.frame(x, Uridylation, status, barcode, position)
      return(dfPercUridylMultiple)
    })
    
    plotUridylationPercentageMultipleReact = reactive({           #Plot pourcentage Uridylation
      plotUridylationPercentageMultiple = ggplot(data = dfPercUridylMultipleReact()) +
        geom_bar(aes(x = x, y = Uridylation, color = as.character(barcode), fill = status), stat = "identity") + 
        facet_wrap(vars(barcode)) + 
        scale_color_hue(name = "", label = paste(listBarcodeMultipleReact())) +
        xlab("Gene (split in different barcodes)") + ylab("Percentage of presence (%)") +
        scale_fill_manual(name = "Additional Tail", values = c("#f1e3ca", "#87E990")) +
        ylim(0, 50) +
        theme(legend.position="bottom", legend.box = "vertical") +
        geom_text(x = 1, y = dfPercUridylMultipleReact()$position, label = paste(round(dfPercUridylMultipleReact()$Uridylation, 2), "%"), size = (5 * (dfPercUridylMultipleReact()$Uridylation/(dfPercUridylMultipleReact()$Uridylation - 0.005))), vjust = 1.05) #le calcul pour la taille permet d'avoo-ir une taille nulle lorsqu'on a 0% mais influe presque pas sur la taille sinon
    })
    
    output$downloadUridylPercMultiplePlot = downloadHandler(        #Télechargement du plot en image
      filename = function(){
        paste0("plotPercentageUridylation-",input$runMultiple, "-", ".png")
      },
      content = function(file){
        ggsave(file, plotUridylationPercentageMultipleReact(), width = 16)
      })
    
    #===================== Plot Nombre U ============================
    
    plotUridylationNumberMultipleReact = reactive({               #Plot nombre de U ajouté 
      plotUridylationNumberMultiple = ggplot(data = dfUridylMultipleReact(), aes(x = x, y = add_tail_T, color = as.character(barcode))) +
        geom_violin(trim = FALSE) +
        facet_wrap(vars(barcode)) +
        geom_boxplot(width = 0.15) +
        ylim(0,10) +
        scale_color_hue(name = "", label = paste(listBarcodeMultipleReact())) +
        theme(legend.position="bottom", legend.box = "vertical") +
        xlab("Gene (split in different barcodes)") + ylab("Number of U")
    })
    
    output$downloadUridylNumberMultiplePlot = downloadHandler(        #Télechargement du plot en image
      filename = function(){
        paste0("plotNumberUridyl-",input$runMultiple, "-", ".png")
      },
      content = function(file){
        ggsave(file, plotUridylationNumberMultipleReact(), width = 16)
      })
    
    # =========================== Plot Distribution taille queue polyA ==============================
    
    
    plotPolyaDistributionMultipleReact = reactive({                 #Plot taille de la queue polyA (densité)
      plotPolyaDistributionMultiple = ggplot(dfUridylMultipleReact(), aes(x = polya_length, group = barcode, color = as.character(barcode))) + geom_density() + xlim(10, 200) +
        scale_color_hue(name = "", labels = paste(listBarcodeMultipleReact())) +
        xlab("PolyA tail length") + ylab("Density of distribution") +
        theme(legend.title = element_text(size = 14), legend.text = element_text(size = 10), legend.position = "bottom")
    })
    
    plotPolyaDistributionMultiple2React = reactive({                #Plot taille de la queue polyA (cumulatif)
      plotPolyaDistributionSingle = ggplot(dfUridylMultipleReact(), aes(polya_length, color = as.character(barcode))) + stat_ecdf(geom = "line") +
        scale_color_hue(name = "", labels = paste(listBarcodeMultipleReact())) +
        theme(legend.title = element_text(size = 14), legend.text = element_text(size = 10), legend.position = "bottom") +
        xlab("PolyA tail length") + ylab("Cumulative proportion") 
    })
    
    plotPolyaDistributionMultiple3React = reactive({                #Plot taille de la queue polyA (boxplot)
      plotPolyaDistributionSingle = ggplot(data = dfUridylMultipleReact(), aes(x = x, y = polya_length, color = as.character(barcode))) +
        geom_violin(trim = FALSE) +
        facet_wrap(vars(barcode)) +
        geom_boxplot(width = 0.15) +
        ylim(0,200) +
        scale_color_hue(name = "", labels = paste(listBarcodeMultipleReact())) +
        theme(legend.position="bottom") +
        xlab("Gene (split in different barcodes)") + ylab("Poly(A) tail length")
    })
    
    output$downloadPolyaDistribMultiplePlot = downloadHandler(           #Télechargement du plot en image, le plot téléchargé dépend de celui affiché donc de la valeur de l'input
      filename = function(){
        paste0("plotDistributionPolya-",input$runMultiple, "-", ".png")
      },
      content = function(file){
        if (input$plotTypeMultiple == "Density plot"){
          ggsave(file, plotPolyaDistributionMultipleReact(), width = 12)}
        else if (input$plotTypeMultiple == "Cumulativ plot"){
          ggsave(file, plotPolyaDistributionMultiple2React(), width = 12)}
        else {
          ggsave(file, plotPolyaDistributionMultiple3React(), width = 12)}
      })
    
    # ===================== Plot Distribution des tailles de queues polyA des mRNA uridylés ==============    
    
    plotPolyaDistributionUridylMultipleReact = reactive({                #Plot taille de la queue polyA (densité)
      polyaDistributionUridylMultiple = ggplot(dfUridylMultipleReact()%>%filter(uridyl == "U"), aes(x = polya_length, group = barcode, color = as.character(barcode))) + geom_density() + xlim(0, 200) +
        scale_color_hue(name = "", labels = paste(listBarcodeMultipleReact())) +
        xlab("PolyA tail length") + ylab("Density of distribution") +
        theme(legend.title = element_text(size = 14), legend.text = element_text(size = 10), legend.position = "bottom")
    })
    
    plotPolyaDistributionUridylMultiple2React = reactive({               #Plot taille de la queue polyA (cumulatif)
      plotPolyaDistributionUridylMultiple = ggplot(dfUridylMultipleReact()%>%filter(uridyl == "U"), aes(polya_length, color = as.character(barcode))) + stat_ecdf(geom = "line") +
        scale_color_hue(name = "", labels = paste(listBarcodeMultipleReact())) +
        theme(legend.title = element_text(size = 14), legend.text = element_text(size = 10), legend.position = "bottom") +
        xlab("PolyA tail length") + ylab("Cumulative proportion") 
    })
    
    plotPolyaDistributionUridylMultiple3React = reactive({               #Plot taille de la queue polyA (boxplot)
      plotPolyaDistributionUridylMultiple = ggplot(data = dfUridylMultipleReact()%>%filter(uridyl == "U"), aes(x = x, y = polya_length, color = as.character(barcode))) +
        geom_violin(trim = FALSE) +
        facet_wrap(vars(barcode)) +
        geom_boxplot(width = 0.15) +
        ylim(0,250) +
        scale_color_hue(name = "", labels = paste(listBarcodeMultipleReact())) +
        theme(legend.position="bottom") +
        xlab("Gene (split in different barcodes)") + ylab("Poly(A) tail length")
    })    
    
    
}

shinyApp(ui, server)










#
#plotUridylationPercentageSingleChangeReact = reactive({
  #plotUridylationPercentageSingle = ggplot(data = dfPercUridylReact()) +
    #geom_bar(aes(x = mRNA, y = Uridylation, fill = status), stat = "identity", position = position_dodge()) + 
    #facet_wrap(vars(barcode)) +
    #xlab("Gene (split in different barcodes)") + ylab("Percentage of presence (%)") +
    #scale_fill_manual(name = "Additional Tail", values = c("#f1e3ca", "#87E990")) +
    #ylim(0, 50) 
#})  
# 