print(sessionInfo())
for (d in c(basilisk::getExternalDir(), tools::R_user_dir("BiocFileCache", "cache"),
            path.expand("~/.cache/huggingface"))) {
    f <- list.files(d, recursive = TRUE, full.names = TRUE, all.files = TRUE)
    cat(d, round(sum(file.size(f), na.rm = TRUE) / 1e9, 2), "GB\n")
}
