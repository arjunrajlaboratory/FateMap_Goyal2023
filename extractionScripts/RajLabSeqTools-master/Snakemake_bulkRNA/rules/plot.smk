rule plotdiffs:
    input: "results/diffexp/{contrast}.diffexp.tsv"
    output:
        volcano="results/diffexp/{contrast}.diffexp.volcano.pdf",
        heatmap="results/diffexp/{contrast}.diffexp.heatmap.pdf"
    threads: 1
    params: 
            numdiffgenes=config["params"]["diffgenes"],
            logsigmax=config["params"]["log-sig"],
            log2foldmin=config["params"]["log2-fold"]
    shell:
        "Rscript scripts/diff-exp-plots.R {input} {params.numdiffgenes} {params.logsigmax} {params.log2foldmin}"
