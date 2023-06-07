# detect system
system = 'sockeye'
node = Sys.info()[["nodename"]]
if (grepl("cedar", node)) {
  system = 'cedar'
} else if (grepl("frontend|worker", node)) {
  system = 'elasti'
} else if (grepl("daint|nid", node)) {
  system = 'cscs'
} else if (Sys.info()["sysname"] == 'Darwin') {
  system = 'local'
}
# set up base directory
if (system == 'cedar') {
  base_dir = paste0("/home/",
                    Sys.info()[["user"]],
                    "/projects/rrg-ljfoster-ab/skinnim/CFdb")
  args$allocation = paste0('rrg-', Sys.info()[["user"]])
} else if (system == 'sockeye') {
  base_dir = "/scratch/st-ljfoster-1/CFdb"
  args$allocation = 'st-ljfoster-1'
} else if (system == 'elasti') {
  base_dir = '/home/ubuntu/projects/rrg-ljfoster-ab/skinnim/CFdb'
  args$allocation = 'root'
} else if (system == 'local') {
  base_dir = "~/git/CFdb-analysis/data"
  args$allocation = NULL
}
