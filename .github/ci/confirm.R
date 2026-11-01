# The Microsoft Store python.exe stubs in WindowsApps ship with every fresh
# Windows install (they open the Store), so they stay: a real user has them.
w <- Sys.which(c("python", "python3", "py", "conda", "git", "uv", "sh", "gcc"))
print(w)
real <- w[c("python", "python3", "py", "conda", "git", "uv")]
stopifnot(!nzchar(real) | grepl("WindowsApps|WINDOW~1", real, ignore.case = TRUE))
