"""Annotate ``altair`` chart with Twitter card and markdown summary.

You can either use function in module or run as command-line tool.

"""


import argparse

from bs4 import BeautifulSoup as bs

import markdown


def annotate_altair_chart(chart_html, annotation_md, twitter_card):
    """
    This function annotates an altair chart with a twitter card and markdown
    description.


    Parameters
    ----------
    chart_html: str
        Path to an HTML file with embeded Vega spec.
    annotation_md: str
        Path to a text file with markdown description to appended to body.
    twitter_card: dict
        Site name, title, description and optionally image for Twitter card.

    Returns
    -------
    str:
        A string of the HTML page formatted to be human readable.
    """

    # Get the main page content
    with open(chart_html, "r") as chart_file:
        page = bs(chart_file, "html.parser")

    # Get the annotation and convert it from markdown to HTML
    with open(annotation_md, "r") as markdow_file:
        annotation = bs(markdown.markdown(markdow_file.read()), "html.parser")

    # Add the annotation to the bottom of the page
    markdown_container = page.new_tag("div", attrs={"id": "markdown"})
    page.body.append(markdown_container)
    separator = page.new_tag("hr")
    markdown_container.append(separator)
    markdown_container.append(annotation)

    # Make and add the twitter card
    if not all(key in twitter_card.keys() for key in ["site", "title", "description"]):
        raise ValueError(
            "Missing required fields for twitter card: site, title, or description"
        )
    summary = page.new_tag("meta", attrs={"name": "twitter:card", "content": "summary"})
    page.head.append(summary)
    for name, content in twitter_card.items():
        card_tag = page.new_tag(
            "meta", attrs={"name": f"twitter:{name}", "content": content}
        )
        page.head.append(card_tag)

    # Add some default styling with bootstrap
    stylesheet = page.new_tag(
        "link",
        attrs={
            "rel": "stylesheet",
            "href": "https://cdn.jsdelivr.net/npm/bootstrap@4.3.1/dist/css/bootstrap.min.css",  # noqa: E501
            "integrity": "sha384-ggOyR0iXCbMQv3Xipma34MD+dH/1fQ784/j6cY/iJTQUOhcWr7x9JvoRxT2MZw1T",  # noqa: E501
            "crossorigin": "anonymous",
        },
    )
    page.head.append(stylesheet)

    # Add some custom styling and margins and overflow
    page.head.style.append(
        "#vis {margin-left: 2.5%; margin-left: 2.5%; width: 95vw; overflow-x: auto;}"
    )
    page.head.style.append(
        "#markdown {margin-left: 2.5%; margin-right: 2.5%; margin-top: 10px; }"
    )
    # Fix the margins and font size for selectors within the vega vis
    page.head.style.append(
        "#vis input, #vis label, #vis span {font-size: 14px; margin: 0px 3px 1px 0px;}"
    )

    return page.prettify()


if __name__ == "__main__":

    # Command line interface
    parser = argparse.ArgumentParser(
        description="Format HTML file containing embeded Vega spec saved with Altair."
    )

    parser.add_argument(
        "--chart",
        type=str,
        required=True,
        help="Path to an HTML file containing a chart saved using Altair.",
    )

    parser.add_argument(
        "--markdown",
        type=str,
        required=True,
        help="Path to a markdown file with text to be included under the plot.",
    )

    parser.add_argument(
        "--site",
        type=str,
        required=True,
        help="URL for the Twitter card.",
    )

    parser.add_argument(
        "--title",
        type=str,
        required=True,
        help="Title of the Twitter card.",
    )

    parser.add_argument(
        "--description",
        type=str,
        required=True,
        help="Description in the Twitter card.",
    )

    parser.add_argument(
        "--image",
        type=str,
        required=False,
        help="Image for Twitter card.",
    )

    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Path to the HTML file to save the formatted chart.",
    )

    args = parser.parse_args()

    # Place the site, title, and description into a dictionary
    twitter_dictionary = {
        "site": args.site,
        "title": args.title,
        "description": args.description,
    }
    if args.image:
        twitter_dictionary["image"] = args.image
    # Get the formated HTML as a string
    annotated_chart = annotate_altair_chart(
        args.chart, args.markdown, twitter_dictionary
    )
    # Write out to a file
    with open(args.output, "w") as outfile:
        outfile.write(annotated_chart)
