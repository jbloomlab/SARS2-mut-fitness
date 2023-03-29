"""Build index for GitHub pages docs."""


import markdown

html_file = snakemake.output.html
plot_annotations = snakemake.params.plot_annotations
mat_trees = snakemake.params.mat_trees
if hasattr(snakemake.params, "current_mat"):
    this_mat = snakemake.params.current_mat
    is_current_mat = True
else:
    this_mat = snakemake.wildcards.mat
    is_current_mat = False

text = [
    f"## {plot_annotations['index_title']}",
    plot_annotations["index_abstract"],
    f"Interactive plots for {this_mat} dataset:",
]

for section, section_title in plot_annotations["sections"].items():
    text.append(f"\n  - {section_title}")
    for plot, plot_d in plot_annotations["plots"].items():
        assert plot_d["section"] in plot_annotations["sections"], plot_d
        if section == plot_d["section"]:
            text.append(
                f"    - [{plot_d['title']}]({plot}.html)"
            )

text.append("\n" + plot_annotations["legend_suffix"])

text.append("\nHere are plots for other datasets:\n")
for mat in mat_trees:
    if is_current_mat:
        text.append(f"  - [{mat} dataset]({mat}/index.html)")
    else:
        text.append(f"  - [{mat} dataset](../{mat}/index.html)")

html = markdown.markdown("\n".join(text))

with open(html_file, "w") as f:
    f.write(html)

