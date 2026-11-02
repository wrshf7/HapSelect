# Tests: simulation progress reporting -----------------------------------------
#
# genomicSimulation is stubbed out, so these pin down only how the simulation loop ticks:
# once per generation, across every rep of both the GA and TS parents.

local_stub_simulation <- function(env = parent.frame()) {
  local_mocked_bindings(
    # Chatty like the real one, which names each file it parses by its path
    load.data            = function(allele.file, map.file, effect.file) {
      cat(sprintf("(Loading %s) 1287 marker(s) with map positions were loaded.\n", map.file))
      list(groupNum = 1L)
    },
    make.random.crosses  = function(...) 1L,
    see.GEBVs            = function(...) c(1, 2),
    break.group.by.GEBV  = function(...) 1L,
    delete.group         = function(...) invisible(NULL),
    see.existing.groups  = function(...) data.frame(Group = integer(0)),
    .package = "genomicSimulation",
    .env = env
  )
}

run_stub_simulation <- function() {
  report_step("local_gebv_simulation")
  run_basic_simulation(as.list(paste0("file", 1:4)), num_gen = 3, num_sim_reps = 2,
                       num_cross_per_gen = 10, num_GA_parents = 2, num_TS_parents = 2,
                       maximize = TRUE)
}

test_that("the simulation ticks once per generation for both sets of parents", {
  local_stub_simulation()
  messages <- character(0)
  events <- with_machine_logging(function() {
    withCallingHandlers(run_stub_simulation(), progression = function(c) {
      if (identical(c$type, "update")) messages <<- c(messages, c$message)
    })
  })
  ticks <- Filter(function(e) identical(e$event, "progress"), events)

  # 2 sets of parents x 2 reps x 3 generations
  expect_length(messages, 12)
  # Count per message rather than order: with_progress() holds conditions back, so they can
  # reach a handler outside it out of order
  expect_equal(as.vector(table(messages)[c("GA parents | rep 1/2", "GA parents | rep 2/2",
                                           "TS parents | rep 1/2", "TS parents | rep 2/2")]),
               c(3L, 3L, 3L, 3L))

  for (tick in ticks) expect_equal(tick$step, "local_gebv_simulation")
  expect_equal(ticks[[length(ticks)]]$current, 12)
  expect_equal(ticks[[length(ticks)]]$total, 12)
})

test_that("the simulation prints what genomicSimulation loaded, labelled rather than by path", {
  local_stub_simulation()
  out <- NULL
  with_machine_logging(function() out <<- capture.output(invisible(run_stub_simulation())))

  # Once per set of parents, with the temporary path swapped for what the file holds,
  # and nothing else: the group deletions stay quiet
  expect_equal(out, c(
    "Loading GA parents into genomicSimulation:",
    "  map:       1287 marker(s) with map positions were loaded.",
    "Loading TS parents into genomicSimulation:",
    "  map:       1287 marker(s) with map positions were loaded."
  ))
})

test_that("the simulation's input files go in the temporary directory, not the working one", {
  one <- data.frame(x = 1)
  files <- write_sim_files(one, one, one, one, temp_string = "test")
  on.exit(delete_sim_files(files), add = TRUE)

  expect_true(all(file.exists(files)))
  expect_equal(unique(normalizePath(dirname(files))), normalizePath(tempdir()))

  delete_sim_files(files)
  expect_false(any(file.exists(files)))
})
