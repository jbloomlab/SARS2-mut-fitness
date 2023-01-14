"""Build index for GitHub pages docs."""


import markdown

html_file = snakemake.output.html
docs_url = snakemake.params.docs_url
plot_annotations = snakemake.params.plot_annotations

text = [
    f"## {plot_annotations['index_title']}",
    plot_annotations["index_abstract"],
    "Interactive plots:",
]

for section, section_title in plot_annotations["sections"].items():
    text.append(f"\n  - {section_title}")
    for plot, plot_d in plot_annotations["plots"].items():
        assert plot_d["section"] in plot_annotations["sections"], plot_d
        if section == plot_d["section"]:
            text.append(
                f"    - [{plot_d['title']}]({docs_url}/{plot}.html)"
            )

text.append("\n" + plot_annotations["legend_suffix"])

html = markdown.markdown("\n".join(text))

with open(html_file, "w") as f:
    f.write(html)

