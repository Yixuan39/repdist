# Stubs every fresh install ships stay, since a real user has them: the
# Microsoft Store python.exe on Windows, /usr/bin/python3 and git on macOS
# (they offer to install the Command Line Tools).
w <- Sys.which(c("python", "python3", "py", "conda", "git", "uv", "sh", "gcc"))
print(w)
real <- w[c("python", "python3", "py", "conda", "git", "uv")]
stub <- grepl("WindowsApps|WINDOW~1", real, ignore.case = TRUE) |
    (Sys.info()[["sysname"]] == "Darwin" & grepl("^/usr/bin/", real))
stopifnot(!nzchar(real) | stub)
