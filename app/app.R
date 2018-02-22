# install.packages("tidyverse")
# install.packages("shiny")
# install.packages("data.table")
# install.packages("seqinr")
# source("https://bioconductor.org/biocLite.R")
# biocLite("ggbio")
# install.packages("shinyBS")
# install.packages("shinyjs")
# install.packages("plotly")
#install.packages("ggiraph")
library(shiny)
library(data.table)
library(plyr)
library(dplyr)
library(seqinr)
library(ggbio)
library(GenomicRanges)
library(plotly)
library(shinyBS)
library(shinyjs)



LangDiamondFiles <- c(paste0("../Similarity Analysis/LangDiamondResults/",c("ERR126028.Allignment.txt","ERR126028_1.Allignment.txt","ERR126028_2.Allignment.txt","ERR126029.Allignment.txt","ERR126029_1.Allignment.txt","ERR126029_2.Allignment.txt")))
SchwarzDiamondFiles <- c(paste0("../Similarity Analysis/SchwarzDiamondResults/",c("SRR935429.Allignment.txt","SRR935429_1.Allignment.txt","SRR935429_2.Allignment.txt")))
LangDiamondFilesBUSCOHMM <- c(paste0("../Similarity Analysis/LangDiamondResultsBuscoHMM/",c("ERR126028.Allignment.txt","ERR126028_1.Allignment.txt","ERR126028_2.Allignment.txt","ERR126029.Allignment.txt","ERR126029_1.Allignment.txt","ERR126029_2.Allignment.txt")))
SchwarzDiamondFilesBUSCOHMM <- c(paste0("../Similarity Analysis/SchwarzDiamondResultsBuscoHMM/",c("SRR935429.Allignment.txt","SRR935429_1.Allignment.txt","SRR935429_2.Allignment.txt")))
LangDiamondFilesMoreSensative <-c(paste0("../Similarity Analysis/LangDiamondResultsMoreSensative/",c("ERR126028.Allignment.txt","ERR126028_1.Allignment.txt","ERR126028_2.Allignment.txt","ERR126029.Allignment.txt","ERR126029_1.Allignment.txt","ERR126029_2.Allignment.txt")))
LangMMseqs2Files <- c(paste0("../Similarity Analysis/LangMMseqs2Results/",c("ERR126028.Allignment.txt","ERR126028_1.Allignment.txt","ERR126028_2.Allignment.txt","ERR126029.Allignment.txt","ERR126029_1.Allignment.txt","ERR126029_2.Allignment.txt")))
  

proteinFileCelegans <-  "../Similarity Analysis/BuscoCelegans.fasta"
proteinFileBusco <-  "../Similarity Analysis/BuscoProteinSet.fa"

Doyle <-  unique(read.delim(file = "../Similarity Analysis/DoyleGenome.tsv",comment.char = "#",header = FALSE,stringsAsFactors = FALSE)[,c(1,2)])
Lang <- unique(read.delim(file = "../Similarity Analysis/LangGenome.tsv",comment.char = "#",header = FALSE,stringsAsFactors = FALSE)[,c(1,2)])
Schwarz <- unique(read.delim(file = "../Similarity Analysis/SchwarzGenome.tsv",comment.char = "#",header = FALSE,stringsAsFactors = FALSE)[,c(1,2)])
colnames(Doyle) <- c("BuscoId","Status")
colnames(Lang) <- c("BuscoId","Status")
colnames(Schwarz) <- c("BuscoId","Status")
Doyle$Status[Doyle$Status=="Complete"] = "Singleton"
Lang$Status[Lang$Status=="Complete"] = "Singleton"
Schwarz$Status[Schwarz$Status=="Complete"] = "Singleton"



diamondFiles <- LangMMseqs2Files
myProteins <- intersect(Doyle$BuscoId[Doyle$Status == "Singleton"],Lang$BuscoId[Lang$Status == "Singleton"])
#myProteins <- Lang$BuscoId[Lang$Status == "Singleton"]
proteinFile <- proteinFileCelegans


#load("LangCelegansBuscoGranges.Rdata")
"LangCelegansBuscoGranges.Rdata"
"SchwarzCelegansBuscoGranges.Rdata"
"SchwarzBuscoGranges.Rdata"
"LangBuscoGranges.Rdata"



readDiamondFiles <- function(diamondFiles,proteinFile,bestAllignment = TRUE){
  
  allDiamondGranges = GRanges()
  for (diamondFile in diamondFiles){
    diamond <- fread(diamondFile,header = FALSE,stringsAsFactors = FALSE,data.table = FALSE)
    colnames(diamond) <- c("Query","Reference","%id","length","mistmatches","gaps","q.start","q.end","s.start","s.end","Diamondevalue","DiamondScore")
    
    if(bestAllignment){
      diamond <-   data.frame(diamond %>%
                                group_by(Query) %>%
                                filter(DiamondScore == max(DiamondScore)))
    }
    diamondGRanges <- makeGRangesFromDataFrame(diamond,seqnames.field= "Reference",start.field = "s.start",end.field ="s.end")
    mcols(diamondGRanges)$Score <- diamond$DiamondScore
    mcols(diamondGRanges)$Query <- diamond$Query
    mcols(diamondGRanges)$file <- strsplit(basename(diamondFile),".",TRUE)[[1]][1]
    allDiamondGranges <-  c(allDiamondGranges,diamondGRanges)
  }
  
  #Correcting sequence lengths 
  if(!is.null(proteinFile)){
    protein <-  read.fasta(proteinFile)
    seqlengths(allDiamondGranges) <-  getLength(protein[       seqnames(seqinfo(allDiamondGranges))        ])
  }
  return(allDiamondGranges)
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

calculateTotalCoverage <- function(coverage){
  totalCoverage <- coverage[,3]
  totalCoverage <- totalCoverage[(totalCoverage > 10 & totalCoverage<=200)]   #Limited the data so it is easier to see the main parts of it. We can also limit the data later on as well. 
  totalCoverage <- as.data.frame(totalCoverage)
  return(totalCoverage)
}


# Define UI for application that draws a histogram
ui <- fluidPage(
    useShinyjs(),
   # Application title
   titlePanel("Sparkle: Making genomes ShinyR"),
   
   #h3("Determine Fold Coverage Threshold"),
   
   
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         sliderInput("Coverage",
                     "Coverage Threshold",
                     min = 10,
                     max = 200,
                     step = 1,
                     value = 0),
         p("Note: A threshold value of '0' does not produce any graphs"),
         actionButton("ResetThresholdCoverage","Reset")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("ThreshHoldCoverageHist")
      )
   ),
   
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
     sidebarPanel(
       
       selectInput("Proteins", "Chosen Proteins", c(), multiple=TRUE, selectize=TRUE),
       sliderInput("AAThreshold","Percent amino acid above threshold",min=0,max=100,value = c(80,100)),
       actionButton("ProtienComplete","Start Analysis"),
       actionButton("myProtein","Select My Proteins")
     ),
     
     # Show a plot of the generated distribution
     mainPanel(
       plotlyOutput("ViolinPlot",width="900",height="500")
       #ggiraphOutput("ViolinPlot",width = "100%", height = "500px")
     )
     

   ),
   

   
   sidebarLayout(
     sidebarPanel(
       selectInput('Protein', 'Protein', c(""), selectize=TRUE),
       p(id="score","Percent amino acid above threshold:"),
       checkboxInput("ShowReads","Show Reads",TRUE),
       checkboxInput("ShowCoverage","Show Coverage",TRUE),
       hr(),
       checkboxInput("BitScoreThreshold","BitScore Threshold",FALSE),
       conditionalPanel("input.BitScoreThreshold == true",checkboxInput("IncludeBadReads","Include Bad Reads",TRUE)),
       conditionalPanel("input.BitScoreThreshold == true",sliderInput("BitScoreThresholdValue","Value",0,100,50)),
       hr(),
       actionButton("Previous","Previous"),
       bsButton("Keep","Keep"),
       bsButton("Discard","Discard"),
       actionButton("Next","Next"),
       br(),
       downloadButton("DownloadResults","Download Results")
     ),
     mainPanel(
       plotlyOutput("ProtienReadGraph",width="900",height="500")
     )
     
   ),
   dataTableOutput("ProteinTable")
   
   
)

# Define server logic required to draw a histogram
server <- function(input, output,session) {
  
  
  #Data 
  
  diamondGRanges <- readDiamondFiles(diamondFiles = diamondFiles,proteinFile = proteinFile)
  save(diamondGRanges,file = "Granges.Rdata")
  splitDiamondGRanges <- split(diamondGRanges,seqnames(diamondGRanges))
  
  coverage <- as.data.frame(coverage(diamondGRanges))
  #Calculating Total Coverage
  totalCoverage <- calculateTotalCoverage(coverage)
  #Make coverage equal to the mode 
  updateSliderInput(session,"Coverage",value =getmode(totalCoverage[,1]) )


  output$ThreshHoldCoverageHist <- renderPlot({
    if(input$Coverage == 0){
      return(NULL)
    }
     ggplot(totalCoverage)+geom_histogram(binwidth = 1,aes(x=totalCoverage,fill= totalCoverage == input$Coverage),show.legend = FALSE) + scale_fill_manual(values = c("FALSE"="black","TRUE"="blue"))+ylab("Number of amino acids")+xlab("Depth of Coverage")
   })
  
  
  
  #Calculaing % of aa above threshold coverage 
  pertcentAAGreaterThanThreshold <- reactive({
    pertcentAAGreaterThanThreshold <-   coverage[,c(2,3)] %>%
      group_by(group_name) %>%
      filter(value >= input$Coverage)%>%
      tally()

    pertcentAAGreaterThanThreshold$Length <- as.vector(seqlengths(diamondGRanges)[pertcentAAGreaterThanThreshold$group_name])
    pertcentAAGreaterThanThreshold <- data.frame(mutate(pertcentAAGreaterThanThreshold,Percentcov=n/Length*100))
    pertcentAAGreaterThanThreshold <- left_join(data.frame(group_name = names(read.fasta(proteinFile))),pertcentAAGreaterThanThreshold)
    pertcentAAGreaterThanThreshold$Percentcov[is.na(pertcentAAGreaterThanThreshold$Percentcov)] = 0
    updateSelectInput(session,"Proteins",choices = pertcentAAGreaterThanThreshold$group_name)
    return(pertcentAAGreaterThanThreshold)
  })
  

  
  ViolinPlot <- reactive({
    plot <- ggplot(pertcentAAGreaterThanThreshold(),aes(x=1,y=Percentcov))+geom_violin(alpha=0.2)
    plot <- plot +coord_flip()
  })
  



  
  output$ViolinPlot <- renderPlotly({
    if(input$Coverage == 0){
      return(NULL)
    }
    set.seed(1)

    plot <- ViolinPlot()
    
    plot <- plot +geom_jitter(shape=16,size=0.4,aes(color=group_name%in%input$Proteins,alpha=ifelse(group_name%in%input$Proteins,100,20),key=group_name,text=group_name))+ scale_colour_manual(values = c("FALSE"="black","TRUE"="blue"))
    plot <- plot +geom_hline(yintercept = isolate(input$AAThreshold),color="black")+ylab(paste0("% amino acid that have coverage above ",input$Coverage))
    ggplotly(plot,width="900",height="500",tooltip=c("text","y"),source="ViolinPlot")%>%layout(showlegend = FALSE)

  })
  
  
  observeEvent(input$AAThreshold+input$Coverage,{
    if(input$Coverage ==0){
      return()
    }
    updateSelectInput(session,"Proteins",selected = c(pertcentAAGreaterThanThreshold()[pertcentAAGreaterThanThreshold()$Percentcov>=input$AAThreshold[1] & pertcentAAGreaterThanThreshold()$Percentcov<=input$AAThreshold[2],1]))
  })
  
  #Select / Click event on the Violoin plot 
  # observeEvent(event_data("plotly_selected",source="ViolinPlot"),{
  #   d <- event_data("plotly_selected",source="ViolinPlot")
  #   if (!is.null(d)){
  #     updateSelectInput(session,"Proteins",selected = c(input$Proteins,d$key))
  #   }
  # })
  
  observeEvent(event_data("plotly_click",source="ViolinPlot"),{
    d <- event_data("plotly_click",source="ViolinPlot")
    if (!is.null(d)){
      updateSelectInput(session,"Proteins",selected = c(input$Proteins,d$key))
    }
  })
  
  # observeEvent(input$ProtienComplete,{
  #   output$ProteinTable <- renderDataTable({
  #     removeUI(selector = "div[id^=Graph]",multiple = TRUE,immediate = TRUE)
  # 
  #     proteinList <- pertcentAAGreaterThanThreshold()
  #     proteinList <- isolate(proteinList[proteinList[,1]%in%input$Proteins,c(1,4)])
  #     colnames(proteinList) <- c("Protein","Percent AA above threshold")
  #     proteinList <- proteinList[order(-proteinList[,2]),]
  #     
  #     sapply(proteinList[,1],graph,splitDiamondGRanges=splitDiamondGRanges,input=input,output=output,session=session)
  #     return(proteinList)
  #     
  #   },options = list(paging = FALSE,searching = FALSE))
  # })
  
  proteinList <- function(){
    proteinList <- pertcentAAGreaterThanThreshold()
    proteinList <- isolate(proteinList[proteinList[,1]%in%input$Proteins,c(1,4)])
    colnames(proteinList) <- c("Protein","Percent AA above threshold")
    proteinList <- proteinList[order(-proteinList[,2]),]
    proteinList$Status <- rep(NA,length(proteinList[,1]))
    return(proteinList)
  }
    
  
  ananlysisProtienList <- reactiveVal()
  
  observeEvent(input$ProtienComplete,{
    
    ananlysisProtienList(proteinList())

    output$ProteinTable <- renderDataTable({
      ananlysisProtienList()
    },options = list(paging = FALSE,searching = FALSE))
    
    updateSelectInput(session,"Protein",choices = c(ananlysisProtienList()[,1]))

  
  })
  
  
  observeEvent(input$Protein,{
    if(input$Protein==""){
      return()
    }
    html("score",paste0("Percent amino acid above threshold: ","<b>",ananlysisProtienList()[ananlysisProtienList()$Protein == input$Protein,2],"</b>"))
    updateButton(session,"Keep",style="default")
    updateButton(session,"Discard",style="default")
    if(is.na(ananlysisProtienList()[ananlysisProtienList()$Protein == input$Protein,3])){
      return()
    }else if(ananlysisProtienList()[ananlysisProtienList()$Protein == input$Protein,3] == "Keep"){
      updateButton(session,"Keep",style="success")
    }else if(ananlysisProtienList()[ananlysisProtienList()$Protein == input$Protein,3]== "Discard"){
      updateButton(session,"Discard",style="danger")
    }

    
  })
  
  observeEvent(input$Keep,{
    if(input$Protein==""){
      return()
    }
    tempProteinList <- ananlysisProtienList()
    tempProteinList[tempProteinList$Protein == input$Protein,3] <- "Keep"
    ananlysisProtienList(tempProteinList)
    
    if(which(ananlysisProtienList()$Protein == input$Protein)+1 > length(ananlysisProtienList()[,1])){
      updateButton(session,"Discard",style="default")
      updateButton(session,"Keep",style="success")
      return()
    }
    updateSelectInput(session,"Protein",selected = ananlysisProtienList()[which(ananlysisProtienList()$Protein == input$Protein)+1,1])
    
  })
  
  observeEvent(input$Next,{
    if(input$Protein==""){
      return()
    }
    
    if(which(ananlysisProtienList()$Protein == input$Protein)+1 > length(ananlysisProtienList()[,1])){
      return()
    }
    updateSelectInput(session,"Protein",selected = ananlysisProtienList()[which(ananlysisProtienList()$Protein == input$Protein)+1,1])
    
  })
  
  observeEvent(input$Previous,{
    if(input$Protein==""){
      return()
    }
    
    if(which(ananlysisProtienList()$Protein == input$Protein)-1 <= 0 ){
      return()
    }
    updateSelectInput(session,"Protein",selected = ananlysisProtienList()[which(ananlysisProtienList()$Protein == input$Protein)-1,1])
    
  })
  
  
  
  observeEvent(input$Discard,{
    if(input$Protein==""){
      return()
    }
    tempProteinList <- ananlysisProtienList()
    tempProteinList[tempProteinList$Protein == input$Protein,3] <- "Discard"
    ananlysisProtienList(tempProteinList)
    if(which(ananlysisProtienList()$Protein == input$Protein)+1 > length(ananlysisProtienList()[,1])){
      updateButton(session,"Keep",style="default")
      updateButton(session,"Discard",style="danger")
      return()
    }
    updateSelectInput(session,"Protein",selected = ananlysisProtienList()[which(ananlysisProtienList()$Protein == input$Protein)+1,1])
  })
  
  
  output$ProtienReadGraph <- renderPlotly({
    if(input$Protein==""){
      return(NULL)
    }

    protienGRanges <- splitDiamondGRanges[[input$Protein]]
    #Adding Color
    BitScoreThreshold = input$BitScoreThresholdValue
    protienGRanges$color <- "blue"
    proteinLength <- length(as.data.frame(coverage(protienGRanges)[input$Protein])[,3])
    
    if(input$BitScoreThreshold){
      protienGRanges$color <- "grey"
      protienGRanges$color[(protienGRanges$Score>=BitScoreThreshold)] <- "blue"
      if(!input$IncludeBadReads){
        protienGRanges <- protienGRanges[protienGRanges$color=="blue"]
      }
    }
    
    readplot <- ggplot()
    if(input$ShowReads){
      readplot <-  readplot +stat_stepping(protienGRanges,aes(color=color,text=Score),geom = "segment",show.legend=FALSE,ylab="Read Depth",xlab="ProtienPosition")+scale_colour_manual(values = c("blue"="blue","grey"="grey")) #Stacked Reads
    }
    if(input$ShowCoverage){
      readplot <- readplot+stat_coverage(protienGRanges,geom="line",color="green",size=2)  #Coverage Alignment
    }
    readplot <- readplot + geom_hline(yintercept = input$Coverage,color="black") + scale_y_continuous()+ coord_cartesian(xlim=c(0,proteinLength))
    ggplotly(readplot,width="900",height="500")%>%layout(showlegend = FALSE)
    
  })
  

  
  output$DownloadResults <- downloadHandler(
    filename="Results.csv",content=function(file) {write.csv(ananlysisProtienList()[ananlysisProtienList()$Status == "Keep",c(1,2)],file)}
  )
  
  #Reset Threshold Coverage
  observeEvent(input$ResetThresholdCoverage,{
    updateSliderInput(session,"Coverage",value =getmode(totalCoverage[,1]))
  })
  
  #Clear Proteins
  # observeEvent(input$ClearProteins,{
  #   updateSelectInput(session,"Proteins",selected = c(""))
  # })
  
  observeEvent(input$myProtein,{
    updateSelectInput(session,"Proteins",selected = myProteins )
  })

  
}

# Run the application 
shinyApp(ui = ui, server = server)

