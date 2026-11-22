source("run_haplotype.R", local = TRUE)

options(shiny.maxRequestSize = max(100 * 1024^2, getOption("shiny.maxRequestSize", 0)))

analysis_defaults <- as.list(formals(run_haplotype))[
  c("max_depth", "split_rule", "min_maf", "fallback_rank",
    "max_missing", "max_heterozygous", "report_max_heterozygous")
]

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
    .haplotype-settings { border-top: 1px solid #ded5e9; margin-top: 20px; padding: 16px 0; }
    .haplotype-settings h2 { margin-top: 0; }
    .haplotype-settings summary { color: #461d7c; cursor: pointer; font-weight: 600; }
    .haplotype-settings details .row { margin-top: 12px; }
    .haplotype-settings .action-button { margin-top: 14px; }
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
      shiny::tags$section(
        class = "haplotype-settings",
        shiny::tags$h2(class = "h4", "Analysis settings"),
        shiny::fluidRow(
          shiny::column(6, shiny::numericInput(
            "max_depth", "Maximum split levels", analysis_defaults$max_depth, min = 1, step = 1
          )),
          shiny::column(6, shiny::selectInput(
            "split_rule", "SNP selection rule",
            choices = c(
              "1. Highest MAF" = "highest_maf",
              "2a. Minimum MAF" = "min_maf",
              "2b. Minimum MAF with fallback" = "min_maf_fallback"
            ),
            selected = analysis_defaults$split_rule
          ))
        ),
        shiny::conditionalPanel(
          "input.split_rule === 'min_maf' || input.split_rule === 'min_maf_fallback'",
          shiny::fluidRow(
            shiny::column(6, shiny::numericInput(
              "min_maf", "Minimum MAF (%)", 100 * analysis_defaults$min_maf, min = 0, max = 50, step = 1
            )),
            shiny::column(6, shiny::conditionalPanel(
              "input.split_rule === 'min_maf_fallback'",
              shiny::numericInput(
                "fallback_rank", "Top MAF markers to consider", analysis_defaults$fallback_rank, min = 1, step = 1
              )
            ))
          ),
          shiny::helpText(paste(
            "Among SNPs meeting the cutoff, the rule uses the smallest count of known homozygous calls.",
            "Fallback considers the top MAF markers when none meets the cutoff."
          ))
        ),
        shiny::tags$details(
          shiny::tags$summary("Marker quality filters"),
          shiny::fluidRow(
            shiny::column(4, shiny::numericInput(
              "max_missing", "Maximum missing calls (%)", 100 * analysis_defaults$max_missing,
              min = 0, max = 100, step = 1
            )),
            shiny::column(4, shiny::numericInput(
              "max_heterozygous", "Maximum unresolved or heterozygous calls for splitting (%)",
              100 * analysis_defaults$max_heterozygous, min = 0, max = 100, step = 1
            )),
            shiny::column(4, shiny::numericInput(
              "report_max_heterozygous", "Maximum unresolved or heterozygous calls in filtered report (%)",
              100 * analysis_defaults$report_max_heterozygous, min = 0, max = 100, step = 1
            ))
          )
        ),
        shiny::actionButton("run_analysis", "Run Analysis", class = "btn-primary"),
        shiny::helpText("Uploads run with the current settings. Click Run Analysis after making changes.")
      ),
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

  analysis_settings <- shiny::reactive({
    number <- function(key, divisor = 1) {
      value <- input[[key]]
      if (is.null(value)) as.numeric(analysis_defaults[[key]]) else as.numeric(value) / divisor
    }
    rule <- if (is.null(input$split_rule)) analysis_defaults$split_rule else input$split_rule
    list(
      max_depth = number("max_depth"),
      split_rule = rule,
      min_maf = if (rule == "highest_maf") analysis_defaults$min_maf else number("min_maf", 100),
      fallback_rank = if (rule == "min_maf_fallback") {
        number("fallback_rank")
      } else {
        analysis_defaults$fallback_rank
      },
      max_missing = number("max_missing", 100),
      max_heterozygous = number("max_heterozygous", 100),
      report_max_heterozygous = number("report_max_heterozygous", 100)
    )
  })

  shiny::observeEvent(list(upload(), selected_sheet(), input$run_analysis), {
    info <- upload()
    if (is.null(info$path)) return()

    sheet <- selected_sheet()
    settings <- analysis_settings()
    key <- list(path = info$path, sheet = sheet, settings = settings)
    if (identical(report()$key, key) && is.null(report()$error)) return()
    clear_report()

    directory <- tempfile("report-", tmpdir = session_dir)
    processed <- shiny::withProgress(message = "Characterizing haplotypes", value = 0.2, {
      tryCatch({
        paths <- do.call(run_haplotype, c(list(
          input_csv = info$path, output_dir = directory, sheet = sheet
        ), settings))
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

  settings_pending <- shiny::reactive({
    result <- report()
    !is.null(result) && !identical(result$key$settings, analysis_settings())
  })

  current_report <- shiny::reactive({
    if (!settings_pending()) report()
  })

  output$run_status <- shiny::renderUI({
    info <- upload()
    if (settings_pending()) {
      return(shiny::tags$p(role = "status", "Settings changed. Click Run Analysis to apply them."))
    }
    result <- current_report()
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
      shiny::helpText(source),
      if (groups == 1L) shiny::helpText("No groups split under these settings.")
    )
  })

  output$downloads <- shiny::renderUI({
    if (is.null(current_report()$paths)) return(NULL)
    shiny::tags$div(
      class = "haplotype-downloads",
      shiny::downloadButton("download_workbook", "Download Excel report", class = "btn-primary"),
      shiny::downloadButton("download_csvs", "Download CSVs")
    )
  })

  output$group_summary <- shiny::renderUI({
    if (is.null(current_report()$paths)) return(NULL)
    shiny::tags$div(
      class = "haplotype-summary",
      shiny::tags$h2(class = "h4", "Haplotype groups"),
      shiny::tableOutput("groups")
    )
  })

  output$groups <- shiny::renderTable(current_report()$groups, rownames = FALSE)

  output$download_workbook <- shiny::downloadHandler(
    filename = function() "haplotype_report.xlsx",
    content = function(file) {
      path <- shiny::req(current_report()$paths[["workbook"]])
      if (!file.copy(path, file, overwrite = TRUE)) stop("The workbook could not be downloaded.")
    }
  )

  output$download_csvs <- shiny::downloadHandler(
    filename = function() "haplotype_csvs.zip",
    content = function(file) {
      result <- current_report()
      shiny::req(result$paths)
      csvs <- result$paths[setdiff(names(result$paths), "workbook")]
      zip::zipr(file, basename(csvs), root = result$directory)
    }
  )
}

shiny::shinyApp(ui, server)
