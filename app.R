library(shiny)
library(shinydashboard)
library(tidyverse)
library(here)
library(DT)
library(plotly)


# Download necessary data -------------------------------------------------
if(!file.exists(here("www", "clover_data.rds"))) {
    viruses <- read_csv(gzcon(url("https://github.com/viralemergence/clover/raw/main/clover/clover_1.0_allpathogens/CLOVER_1.0_Viruses_AssociationsFlatFile.csv")), show_col_types = FALSE)
    
    arenaviridae_hantaviridae <- viruses %>%
        filter(str_detect(PathogenFamily, "arenaviridae|hantaviridae")) %>%
        filter(str_detect(HostOrder, "rodentia")) %>%
        select(AssocID, Host, HostGenus, Pathogen, DetectionMethod, DatabaseVersion, DatabaseDOI, HostTaxID, PathogenTaxID) %>%
        mutate(Reference = paste0("<a href='", DatabaseDOI, "'>", DatabaseVersion, "</a>"),
               Host = str_to_sentence(Host),
               HostGenus = str_to_sentence(HostGenus),
               Pathogen = str_to_sentence(Pathogen)) %>%
        group_by(Host, Pathogen, DetectionMethod) %>%
        arrange(AssocID) %>%
        slice(1) %>%
        select("Association ID" = AssocID,
               "Host species" = Host,
               "Host genus" = HostGenus,
               Pathogen,
               "Detection method" = DetectionMethod,
               "Host ID" = HostTaxID,
               "Pathogen ID" = PathogenTaxID,
               "Database reference" = DatabaseVersion,
               "Database DOI" = DatabaseDOI,
               Reference)
    
    write_rds(arenaviridae_hantaviridae, here("www", "clover_data.rds"))
    
}

clover_data <- read_rds(here("www", "clover_data.rds"))

# UI ----------------------------------------------------------------------

# Define UI for application that draws a histogram
ui <- dashboardPage(
    
    skin = "yellow",
    
    dashboardHeader(title = "Arenaviruses and Hantaviruses of rodents",
                    titleWidth = 450),
    
    dashboardSidebar(
        width = 350,
        sidebarMenu(
            menuItem("Homepage", tabName = "homepage", icon = icon("house")),
            menuItem("Known Arenaviridae/Hantaviridae of rodents", tabName = "known_pathogens", icon = icon("viruses")),
            menuItem("Rodent hosts", tabName = "rodent_range", icon = icon("mouse"))
        )
    ),
    
    
    ## Homepage ----------------------------------------------------------------
    
    dashboardBody(
        tabItems(
            # Homepage tab content
            tabItem(tabName = "homepage",
                    fluidRow(
                        
                        box(width = 12,
                            h1("Introduction"),
                            p("Rodents are global hosts of zoonotic pathogens and potential hosts of novel pathogens of epidemic potential. Existing efforts to catalogue host-pathogen associations in these species are limited by global datasets which lack temporal and geographic specificity. Current research is hindered by spatial-, host taxa- and temporal biases within these datasets that are challenging to quantify. Here, we produce a database of studies on rodent-pathogen associations on a global scale focussed on two important viral families, ", strong("Arenaviridae"), " and ", strong("Hantaviridae"), "."),
                            p("Arenaviridae and Hantaviridae are two important, globally distributed, rodent-associated viral families in the order ", strong("Bunyavirales"), ". This order contains several known zoonoses. Arenaviridae associated zoonoses include ", em("Lassa mammarenavirus"), " the cause of Lassa fever in West Africa, ", em("Lymphocytic choriomeningitis virus"), " the cause of Lymphocytic choriomeningitis disease. Hantaviridae associated zoonoses include, ", em("Hantaan orthohantavirus"), " and ", em("Sin Nombre orthohantavirus"), " which cause Hantavirus haemorrhagic fever with renal syndrome and Hantavirus pulmonary syndrome respectively."),
                            p("This data contained within this application has been extracted from a systematic review of the peer-reviewed scientific literature, pre-print articles, ecological reports and `grey` literature with data extracted for subsequent analysis. Data has been extracted on rodent sampling, including species sampled, sampling locations and sampling effort, which is matched to pathogen assays. Where possible these assays have been linked to pathogen sequences archived on NCBI GenBank and pathogenesis studies conducted on these viruses."),
                            p("These data are made available in several formats:"),
                            # Numbered bullet points for the data contained within the application from output$pageslist
                            uiOutput("pageslist"),
                            p(HTML("The repository for the data extraction project and the study protocol is available on GitHub within the <a href=`https://github.com/DidDrog11/arenavirus_hantavirus`>arenavirus_hantavirus</a> project."))
                        )
                    )),
            
            
            ## Known pathogens ---------------------------------------------------------
            
            # Known Arenavirus and Hantavirus pathogens of rodents
            tabItem(tabName = "known_pathogens",
                    fluidRow(
                        
                        box(width = 12,
                            h1("Datasource"),
                            p("CLOVER is a harmonised mammal-virus association database constructed to reconcile and aggregate four datasets that catalogue host-virus associations", HTML("<a href=`https://doi.org/10.1093/biosci/biab080`>(Gibb et al. 2021)</a>"), ". The table below shows all of the contained host-pathogen associations for species of Rodentia and Arenaviridae or Hantaviridae.")),
                        
                        box(width = 12,
                            h1("Table of host-pathogen associations for Rodentia and Arenaviridae or Hantaviridae"),
                            p("151 rodent host species of 48 pathogens are contained within the CLOVER database. 24 of these pathogens are Arenaviridae and 24 Hantaviridae."),
                            downloadButton("download_clover", "Download table as .csv"),
                            br(),
                            p("This will download the entire dataset or a subset depending on the results of the search terms in the search box."),
                            # This datatable contains the output of the clover_data rodent-arenaviridae/hantaviridae subset
                            DTOutput("clover_data"))
                    )),
            
            
            ## Rodent ranges -----------------------------------------------------------
            
            # Ranges of rodents known to be hosts of Arenaviridae and Hantaviridae
            tabItem(tabName = "rodent_range",
                    fluidRow(
                        
                        box(width = 12,
                            h1("Geographic ranges of known rodent hosts of Arenaviridae and Hantaviridae"),
                            p("151 rodent species are identified are identified to be hosts of pathogens. Most rodent species are hosts of a single pathogen (70, 46%), with few species host to 5 or more pathogen species (10%)."),
                        ),
                        
                        box(width = 6,
                            plotlyOutput("rodent_hosts"))
                        
                        
                    ))
            
        )
    )
    
)

# Server ------------------------------------------------------------------

# Define server logic required to draw a histogram
server <- function(input, output) {
    

# Homepage ----------------------------------------------------------------

    # Numbered bullet points for the data contained within the application
    output$pageslist <- renderUI({
        tags$ol(
            tags$li("The included studies with citations and links to the original publication from which data were extracted."),
            tags$li("A description of the rodent species sampled for pathogens of these two viral families."),
            tags$li("The location of studies and sampling of rodents."),
            tags$li("The prevalence of pathogens of these viral families among sampled rodent species."),
            tags$li("The location of detections of these viruses within hosts."),
            tags$li("Locations of available viral sequences for detected pathogens.")
        )
    })
    

# Known pathogens ---------------------------------------------------------
    
    # Render a datatable of the clover dataframe that is stored in the www folder and loaded in as clover_data
    output$clover_data <- renderDT(
        datatable(clover_data %>%
                      select("Association ID",
                             "Host species",
                             Pathogen,
                             "Detection method",
                             "Host ID",
                             "Pathogen ID",
                             Reference),
                  escape = FALSE)
    )
    
    # This produces the filtered data based on the search field
    output$filtered_clover <- renderPrint({
        input[["clover_data_rows_all"]]
    })
    
    # This button downloads the filtered data
    output$download_clover <- downloadHandler(filename = "clover_rodent_pathogen.csv",
                                              content = function(file) {
                                                  write.csv(clover_data[input[["clover_data_rows_all"]], ] %>%
                                                                select(-Reference),
                                                            file,
                                                            row.names = FALSE)
                                                  })
    

# Rodent ranges -----------------------------------------------------------
    
    output$rodent_hosts <- renderPlotly({
        
        plot_rodent <- clover_data %>%
            ungroup() %>%
            select(`Host species`, `Host genus`, Pathogen) %>%
            group_by(`Host species`) %>%
            mutate(`N pathogens` = n()) %>%
            group_by(`Host genus`) %>%
            mutate(species_number = match(`Host species`, unique(`Host species`))) %>%
            ungroup() %>%
            mutate(`Host genus` = fct_rev(fct_infreq(`Host genus`)))
        
        
        plot_ly(plot_rodent,
                x = ~`N pathogens`,
                y = ~`Host genus`,
                color = ~factor(species_number),
                stroke = I("black"),
                colors = "YlOrRd",
                type = "bar") %>%
            layout(barmode = "stack")
        
    })


    
}


# Run application ---------------------------------------------------------

# Run the application 
shinyApp(ui = ui, server = server)
