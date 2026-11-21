source("run_haplotype.R", local = TRUE)

options(shiny.maxRequestSize = max(100 * 1024^2, getOption("shiny.maxRequestSize", 0)))

ui <- shiny::fluidPage(
  title = "Haplotype Characterization",
  shiny::tags$head(shiny::tags$style(shiny::HTML("
    body { background: #f5f3f8; color: #231b32; }
    .haplotype-page { max-width: 780px; margin: 56px auto; padding: 0 12px; }
    .haplotype-page h1 { color: #461d7c; font-weight: 700; }
    .haplotype-page .intro { color: #50465d; font-size: 17px; margin-bottom: 24px; }
    .haplotype-panel {
      background: white; border: 1px solid #ded5e9; border-top: 4px solid #461d7c;
      border-radius: 12px; padding: 24px; box-shadow: 0 12px 32px #2815420c;
    }
    .haplotype-panel .btn-primary { background: #461d7c; border-color: #461d7c; }
    .haplotype-panel .btn-primary:hover { background: #35145e; border-color: #35145e; }
    .haplotype-downloads { display: flex; flex-wrap: wrap; gap: 8px; margin: 14px 0; }
    .haplotype-summary { margin-top: 22px; }
    .haplotype-summary table { max-width: 420px; }
    @media (max-width: 600px) {
      .haplotype-page { margin: 26px auto; }
      .haplotype-panel { padding: 16px; }
    }
  "))),
  shiny::tags$main(
    class = "haplotype-page",
    shiny::tags$h1("Haplotype Characterization"),
    shiny::tags$p(class = "intro", "Upload SNP calls to generate the haplotype report."),
    shiny::tags$section(
      class = "haplotype-panel",
      shiny::fileInput("snp_file", "SNP file", accept = c(".csv", ".xlsx"), width = "100%"),
      shiny::helpText("Use a Line or germplasmName column and SNP_chr marker columns."),
      shiny::uiOutput("sheet_picker"),
      shiny::uiOutput("run_status"),
      shiny::uiOutput("downloads"),
      shiny::uiOutput("group_summary")
    )
  )
)

server <- function(input, output, session) {
  session_dir <- tempfile("haplotype-session-")
  if (!dir.create(session_dir)) stop("Could not create a temporary directory for this session.")
  session$onSessionEnded(function() unlink(session_dir, recursive = TRUE))

  upload <- shiny::reactiveVal(NULL)
  report <- shiny::reactiveVal(NULL)

  clear_report <- function() {
    previous <- report()
    if (!is.null(previous$directory)) unlink(previous$directory, recursive = TRUE)
    report(NULL)
  }

  shiny::observeEvent(input$snp_file, {
    previous <- upload()
    if (!is.null(previous$directory)) unlink(previous$directory, recursive = TRUE)
    clear_report()

    file <- input$snp_file
    name <- basename(file$name)
    extension <- tolower(tools::file_ext(name))
    if (!extension %in% c("csv", "xlsx")) {
      upload(list(error = "Upload a CSV or XLSX file."))
      return()
    }

    directory <- tempfile("snp-", tmpdir = session_dir)
    path <- file.path(directory, name)
    info <- tryCatch({
      if (!dir.create(directory) || !file.copy(file$datapath, path)) {
        stop("Could not read the uploaded file.", call. = FALSE)
      }
      sheets <- if (extension == "xlsx") haplotype_input_sheets(path) else character()
      if (extension == "xlsx" && !length(sheets)) {
        stop("The workbook needs a worksheet with a Line or germplasmName column and SNP_chr markers.",
             call. = FALSE)
      }
      list(path = path, directory = directory, name = name, sheets = sheets)
    }, error = function(e) {
      unlink(directory, recursive = TRUE)
      list(error = conditionMessage(e))
    })
    upload(info)
  }, priority = 10)

  selected_sheet <- shiny::reactive({
    sheets <- upload()$sheets
    selected <- intersect(input$sheet, sheets)
    if (length(selected)) selected[1L] else if (length(sheets)) sheets[1L]
  })

  output$sheet_picker <- shiny::renderUI({
    sheets <- upload()$sheets
    if (length(sheets) > 1L) {
      shiny::selectInput(
        "sheet", "SNP worksheet", choices = sheets,
        selected = shiny::isolate(selected_sheet())
      )
    }
  })

  shiny::observeEvent(list(upload(), selected_sheet()), {
    info <- upload()
    if (is.null(info$path)) return()

    sheet <- selected_sheet()
    key <- c(info$path, sheet)
    if (identical(report()$key, key)) return()
    clear_report()

    directory <- tempfile("report-", tmpdir = session_dir)
    processed <- shiny::withProgress(message = "Characterizing haplotypes", value = 0.2, {
      tryCatch({
        paths <- run_haplotype(info$path, directory, sheet = sheet)
        summary <- utils::read.csv(paths[["representatives"]], check.names = FALSE)
        list(
          key = key, directory = directory, paths = paths, sheet = sheet,
          groups = data.frame(Group = summary$Concatenated_Grouping, Lines = summary$Grouping_Count)
        )
      }, error = function(e) {
        unlink(directory, recursive = TRUE)
        message <- sub(info$path, info$name, conditionMessage(e), fixed = TRUE)
        list(key = key, error = message)
      })
    })
    report(processed)
  })

  output$run_status <- shiny::renderUI({
    info <- upload()
    result <- report()
    error <- if (!is.null(info$error)) info$error else result$error
    if (!is.null(error)) {
      return(shiny::tags$p(class = "text-danger", role = "alert", error))
    }
    if (is.null(result$paths)) return(NULL)

    source <- if (is.null(result$sheet)) info$name else paste(info$name, result$sheet, sep = " · ")
    lines <- sum(result$groups$Lines)
    groups <- nrow(result$groups)
    shiny::tagList(
      shiny::tags$p(
        role = "status", shiny::tags$strong("Report ready. "),
        sprintf(
          "%s in %s.", ngettext(lines, paste(lines, "line"), paste(lines, "lines")),
          ngettext(groups, paste(groups, "group"), paste(groups, "groups"))
        )
      ),
      shiny::helpText(source)
    )
  })

  output$downloads <- shiny::renderUI({
    if (is.null(report()$paths)) return(NULL)
    shiny::tags$div(
      class = "haplotype-downloads",
      shiny::downloadButton("download_workbook", "Download Excel report", class = "btn-primary"),
      shiny::downloadButton("download_csvs", "Download CSVs")
    )
  })

  output$group_summary <- shiny::renderUI({
    if (is.null(report()$paths)) return(NULL)
    shiny::tags$div(
      class = "haplotype-summary",
      shiny::tags$h2(class = "h4", "Haplotype groups"),
      shiny::tableOutput("groups")
    )
  })

  output$groups <- shiny::renderTable(report()$groups, rownames = FALSE)

  output$download_workbook <- shiny::downloadHandler(
    filename = function() "haplotype_report.xlsx",
    content = function(file) {
      path <- shiny::req(report()$paths[["workbook"]])
      if (!file.copy(path, file, overwrite = TRUE)) stop("The workbook could not be downloaded.")
    }
  )

  output$download_csvs <- shiny::downloadHandler(
    filename = function() "haplotype_csvs.zip",
    content = function(file) {
      result <- report()
      shiny::req(result$paths)
      csvs <- result$paths[setdiff(names(result$paths), "workbook")]
      zip::zipr(file, basename(csvs), root = result$directory)
    }
  )
}

shiny::shinyApp(ui, server)
