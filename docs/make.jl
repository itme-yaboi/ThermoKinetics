# Build the HTML documentation:  julia --project=docs docs/make.jl   (from the repo root)
using Documenter
using ThermoKinetics

const DOCS = @__DIR__
const PROJ = joinpath(DOCS, "src", "project")
mkpath(PROJ)

# Copy the project records into docs/src/project/ (generated, git-ignored). Markdown links whose
# target is a local file are turned into plain text: the targets live outside the documentation
# tree (or are not published), and Documenter's link checks would fail on them.
function plain_local_links(text)
    # [label](target) with a target that is not http(s), mailto or a same-page anchor
    re = r"(?<!\!)\[([^\]]*)\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)"
    return replace(text, re => function (m)
        mt = match(re, m)
        label, target = mt.captures
        (startswith(target, "http://") || startswith(target, "https://") ||
         startswith(target, "mailto:") || startswith(target, "#")) && return m
        return label
    end)
end

for f in ("DESIGN.md", "STATUS.md", "DECISIONS.md")
    write(joinpath(PROJ, f), plain_local_links(read(joinpath(DOCS, f), String)))
end

makedocs(
    sitename = "ThermoKinetics.jl",
    modules = [ThermoKinetics],
    format = Documenter.HTML(prettyurls = false, size_threshold_warn = 400_000, size_threshold = 1_000_000),
    checkdocs = :exports,
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting_started.md",
        "Architecture and structure" => "architecture.md",
        "API reference" => [
            "Domains" => "api/domains.md",
            "Materials" => "api/materials.md",
            "Nodes" => "api/nodes.md",
            "Conductors and dissipation" => "api/conductors.md",
            "Converters and transformers" => "api/converters.md",
            "Flow" => "api/flow.md",
            "Sources, reactions, nuclear" => "api/sources.md",
            "Switching" => "api/switching.md",
            "System building and running" => "api/system.md",
            "Results and diagnostics" => "api/results.md",
            "Spatial helpers" => "api/spatial.md",
        ],
        "Project records" => [
            "Design" => "project/DESIGN.md",
            "Status" => "project/STATUS.md",
            "Decisions" => "project/DECISIONS.md",
        ],
    ],
)
