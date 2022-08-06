const buildGraphQuilt = (graph, slug) => {
  const graphId = `#allele_graph_${slug}`
  const tableId = `#allele_table_${slug}`
  const assignNumPercentReadsToLinks = link => {
    if(!link.hasOwnProperty("numReads") && !link.hasOwnProperty("percentReads")) {
      return Object.assign(link, {
        numReads: Math.min(
          graph.nodes[link.source].numReads,
          graph.nodes[link.target].numReads,
        ),
        percentReads: Math.min(
          graph.nodes[link.source].percentReads,
          graph.nodes[link.target].percentReads,
        ),
      })
    } else {
      return link
    }
  }
  graph["links"] = graph["links"].map(assignNumPercentReadsToLinks)

  const color = {
    "A": "#7FC97F",
    "T": "#BEAED4",
    "C": "#FDC086",
    "G": "#FFFF99",
    "N": "#C8C8C8",
    "Ø": "#C8C8C8",
    "0": "#C8C8C8",
    "-": "#C1C1C1",
    "highlight": "#E84855"
  }
  // b8336a // maroon X 11
  // e84855 // sizzling red
  // ed1c24 // red pigment
  // 348aa7 // blue munsell
  // 26428B // dark cornflower blue

  const seqLength = graph["alleles"][0].seq.length

  const quiltSvg = d3.select(tableId)
        .append("svg")
        .attr("viewBox", `0 -10 ${(seqLength * 10) + 150} ${graph["alleles"].length * 10}`)
        .style("user-select", "none")


  const alleleGroup = quiltSvg.selectAll(`.${slug}_allele`)
        .data(graph["alleles"], d => d.id)
        .join(
          enter => {
            let root = enter.append("g")
            root.attr("fill", "#999999")
              .attr("stroke", "#999999")
              .attr("stroke-width", "0.5pt")
              .style("font-size", "8pt")
              .style("opacity", 1)
              .attr("transform", (d, i) => `translate(0,${i * 10})`)
              .attr("id", d => `${slug}_allele_${d.id}`)
            root.append("text")
              .attr("x", (seqLength * 10) + 2)
              .attr("dy", "-0.2em")
              .text(d => `${d.name}, ${d.numReads} (${d.percentReads.toFixed(2)}%) Reads`)
            return root
          },
        )

  const alleleSeqGroup = alleleGroup.selectAll(`.${slug}_alleleSeq`)
        .data(d => d.seq, d => d)
        .join(enter => {
          let root = enter.append("g")
          root.append("rect")
            .attr("fill", d => color[d])
            .attr("width", 10)
            .attr("height", 10)
            .attr("x", (d, i) => i * 10)
            .attr("y", -10)
          root.append("text")
            .attr("x", (d, i) => (i * 10) + 5)
            .attr("text-anchor", "middle")
            .attr("dy", "-0.2em")
            .style("font-size", "6pt")
            .style("stroke-width", "0.3pt")
            .text(d => d)
          return root
        })

  // bold face the substitutions
  alleleGroup.filter(d => d.hasOwnProperty("substitutions") && d["substitutions"].length > 0)
    .each(function(d) {
      d3.select(this)
        .selectAll("g")
        .select("text")
        .filter((alleleBase, baseIndex) => d["substitutions"].includes(baseIndex))
        .style("font-weight", "bold")
    })

  // add rectangles around insertions
  alleleGroup.filter(d => d.hasOwnProperty("insertions") && d["insertions"].length > 0)
    .each(function(d) {
      d["insertions"].forEach(e => d3.select(`${tableId} > svg`)
                              .insert("rect")
                              .attr("width", e[1] * 10)
                              .attr("height", 10)
                              .attr("x", (e[0] + 1) * 10)
                              .attr("y", (d.id - 1)* 10)
                              .style("stroke", color["highlight"])
                              .style("fill", "none")
                             )
    })

  let inactiveAlleles = new Set()

  const hideOrShowAlleleGroup = function(d, el, update = false) {
    let newAlleleOpacity = 1,
        newGraphOpacity = 1
    var alleleElements = d3.select(graphId)
        .selectAll(`.${slug}_allele_${d.id}`)
        .data()
    // check if this allele is currently hidden
    if(!inactiveAlleles.has(d.id)) {
      newAlleleOpacity = 0.3
      newGraphOpacity = 0
      inactiveAlleles.add(d.id)
      alleleElements = alleleElements.filter(e => e.alleleIds.filter(f => inactiveAlleles.has(f)).length == e.alleleIds.length)
    }
    else {
      inactiveAlleles.delete(d.id)
    }
    d3.select(el).style("opacity", newAlleleOpacity)
    alleleElements.forEach(e => {
      if(e.hasOwnProperty("id")) {
        graph.nodes[e.id].show = !e.show
      }
    })
    if(update) {
      updateGraph(svg, d3cola)
    }
  }

  const hoverAlleleGroup = (d, nodeStrokeColor, linkStrokeColor) => {
    if(!inactiveAlleles.has(d.id)) {
      d3.selectAll(d.alleleNodeIds.map(e => `#${slug}_node_${e}`).join(","))
        .select("rect")
        .style("stroke", nodeStrokeColor)
      d3.selectAll(`.${slug}_link.${slug}_allele_${d.id}`)
        .style("stroke", linkStrokeColor)
    }
  }

  alleleGroup
    .filter(d => d.name != "Reference")
    .on("click", function(d) { hideOrShowAlleleGroup(d, this, true) })
    .on("mouseover", d => { hoverAlleleGroup(d, color["highlight"], color["highlight"]) })
    .on("mouseout", d => { hoverAlleleGroup(d, "#ffffff", "#999999") })

  const hideOrShowAllAlleleGroups = () => {
    alleleGroup
      .filter(d => d.name != "Reference")
      .each(function(d) { hideOrShowAlleleGroup(d, this, false) })

    updateGraph(svg, d3cola)
  }

  const width = seqLength * 40,
        height = seqLength * 5

  var svg = d3.select(graphId)
      .append("svg")
      .attr("viewBox", `0 0 ${width} ${height}`)

  // add a background to the entire svg
  svg.append("rect")
    .attr("width", "100%")
    .attr("height", "100%")
    .attr("fill", "white")
    .style("opacity", 0.7)

  const zoom = d3.zoom()
        .scaleExtent([1 / (width / window.innerWidth) * 0.75, 20])
        .translateExtent([[0, 0], [width / 2, height / 8]])
        .on("zoom", () => {
          let transform = d3.event.transform
          transform.x = 0
          transform.y = 0
          svg.attr("transform", transform.toString())
        })
  svg.call(zoom)

  let d3cola = cola.d3adaptor(d3)
      .size([width, height])
      .avoidOverlaps(true)

  const maxFrequency = Math.max(graph["alleles"].map(d => d.percentReads)),
        minFrequency = Math.min(graph["alleles"].map(d => d.percentReads)),
        minWidth = 1,
        widthMultiplier = 1.8

  const frequencyExtent = d3.extent(graph["alleles"].map(d => d.percentReads))
  const opacityScale = d3.scaleSequentialPow()
        .domain(frequencyExtent)
        .interpolator(d3.interpolateNumber(0.4, 0.9))

  const widthScale = d3.scaleLinear()
        .domain(frequencyExtent)
        .range([0.05, 0.15])

  const calculateLinkWidth = (frequency) => ((frequency / 100) * widthMultiplier) + minWidth

  const isVisibleNode = node => node.show

  const isInvisibleNode = node => !isVisibleNode(node)

  const isVisibleLink = link => {
    if(typeof(link.source) == "number" && typeof(link.target) == "number") {
      return graph.nodes[link.source].show && graph.nodes[link.target].show
    } else {
      return graph.nodes[link.source.id].show && graph.nodes[link.target.id].show
    }
  }

  const isInvisibleLink = link => !isVisibleLink(link)

  const updateSetcola = (d3cola) => {
    graph.groups = graph.groups.map((group, group_index) => {
      // set the group index for each node
      group.leaves.forEach(d => {
        if(typeof(d) == "number") {
          graph.nodes[d]["group_index"] = group_index
        } else {
          graph.nodes[d.id]["group_index"] = group_index
        }
      })
      // remove the nodes that aren't visible from this group
      group.leaves = group.leaves.filter(d => {
        if(typeof(d) == "number") {
          return graph.nodes[d].show
        } else {
          return d.show
        }
      })
      return group
    })
    // add the nodes that are visible back to the group
    graph.nodes.filter(d => "group_index" in d && d.show)
      .forEach(d => graph.groups[d.group_index].leaves.push(d.id))
    const setcolaResult = setcola
          .nodes(graph.nodes)
          .links(graph.links.filter(isVisibleLink))
          .groups(graph.groups)
          .constraints(graph.setcolaSpec)
          .layout()

    d3cola
      .nodes(setcolaResult.nodes)
      .links(setcolaResult.links)
      .groups(setcolaResult.groups)
      .constraints(setcolaResult.constraints)
      .avoidOverlaps(true)
      .start(10, 10, 10)

    return setcolaResult
  }

  const collapseNode = (node, update = true) => {
    let lastLink = graph.links.filter(link => {
      if(node.collapsedNodes && link.source.id === node.collapsedNodes.slice(-1)[0]) {
        return true
      } else {
        return false
      }
    })[0]
    if(!node.collapsed) {
      node.collapsedNodes
        .slice(1) // skip the first node of the group
        .forEach(e => {
          graph.nodes[e].show = false
          graph.nodes[e].collapsed = true
        })
      graph.nodes[node.collapsedNodes[0]].collapsed = true
      // add a new edge from the collapsed node to the first uncollapsed node
      if(lastLink) {
        graph.links.push({
          type: lastLink.type,
          alleleIds: lastLink.alleleIds,
          numReads: lastLink.numReads,
          percentReads: lastLink.percentReads,
          source: graph.nodes[node.collapsedNodes[0]],
          target: lastLink.target,
        })
        graph.nodes[node.collapsedNodes[0]].extraLink = {
          source: graph.nodes[node.collapsedNodes[0]].id,
          target: lastLink.target.id
        }
      }
    } else {
      node.collapsedNodes
        .forEach(e => {
          graph.nodes[e].show = true
          graph.nodes[e].collapsed = false
        })
      if(node.extraLink) {
        graph.links = graph.links.filter(d => !(d.source.id === node.extraLink.source && d.target.id === node.extraLink.target))
      }
    }
    if(update) {
      updateGraph(svg, d3cola)
    }
  }

  const collapseAllNodes = () => {
    graph.nodes
      .filter(d => isVisibleNode(d) && d.collapsedNodes && !d.collapsed && (d.id == d.collapsedNodes[0]))
      .forEach(d => collapseNode(d, false))
    updateGraph(svg, d3cola)
  }

  const bounds = (value) => Math.min(Math.max(30, value), Math.max(width, height))

  const updateGraph = (svg, d3cola) => {
    let setcolaResult = updateSetcola(d3cola)

    let group = createGroup(svg, d3cola, setcolaResult)
    let groupLabel = createGroupLabel(svg, setcolaResult)
    let link = createLink(svg, setcolaResult)

    let deletions = createDeletions(svg, setcolaResult)


    if(setcolaResult.nodes.filter(d => d.cleavagePosition).length > 0) {
      createCleavagePosition(svg, setcolaResult)
      createCleavageLabel(svg)
    }

    let nodeGroup = createNodeGroup(svg, setcolaResult)
    nodeGroup.raise()

    nodeGroup.filter(d => d.collapsedNodes)
      .on("click", d => {collapseNode(d, true)})

    d3cola.on("tick", () => {
      link
        .attr("x1", d => bounds(d.source.x))
        .attr("y1", d => bounds(d.source.y))
        .attr("x2", d => bounds(d.target.x))
        .attr("y2", d => bounds(d.target.y));

      deletions
        .attr("d", d => {
          var dx = d.target.x - d.source.x,
              dy = d.target.y - d.source.y,
              dr = Math.sqrt(dx * dx + dy * dy);
          return `M${bounds(d.source.x)},${bounds(d.source.y)}A${dr},${dr} 0 0,1 ${bounds(d.target.x)},${bounds(d.target.y)}`
        })

      nodeGroup.select("rect")
        .attr("x", d => bounds(d.x) - d.width / 2)
        .attr("y", d => bounds(d.y) - nodeHeight / 2)

      let group_padding = nodeHeight * 2
      group
        .attr("x", d => {
          d3.select(`#${slug}_${d.name}_label`)
            .attr("x", bounds(d.bounds.x))
          return bounds(d.bounds.x) - group_padding / 4
        })
        .attr("y", d => {
          d3.select(`#${slug}_${d.name}_label`)
            .attr("y", bounds(d.bounds.y) + d.bounds.height() + group_padding)
          return bounds(d.bounds.y) - group_padding / 4
        })
        .attr("width", d => d.bounds.width() + group_padding / 2)
        .attr("height", d => d.bounds.height() + group_padding / 2)

      nodeGroup.select("text")
        .attr("x", d => bounds(d.x))
        .attr("y", d => bounds(d.y) + nodeHeight / 4)
    });
  }

  const getLinkId = link => {
    if(typeof(link.source) == "number" && typeof(link.target) == "number") {
      return `${slug}_link_${link.source}_${link.target}`
    } else {
      return `${slug}_link_${link.source.id}_${link.target.id}`
    }
  }

  const createGroup = (svg, d3cola, setcolaResult) => svg.selectAll(`.${slug}_group`)
        .data(setcolaResult.groups)
        .join("rect")
        .attr("rx", 8)
        .attr("ry", 8)
        .attr("class", `${slug}_group`)
        .style("fill", "#888888")
        .style("opacity", 0.5)
        .call(d3cola.drag)

  const createGroupLabel = (svg, setcolaResult) => svg.selectAll(`.${slug}_groupLabel`)
        .data(setcolaResult.groups)
        .join("text")
        .attr("id", d => `${slug}_${d.name}_label`)
        .attr("class", `${slug}_groupLabel`)
        .text(d => d.name)

  const createLink = (svg, setcolaResult) => svg.selectAll(`.${slug}_link`)
        .data(setcolaResult.links.filter(d => d.type !== "Deletion"), d => getLinkId(d))
        .join(
          enter => {
            let link = enter.append("line")
                .attr("class", d => {
                  if(d.alleleIds) {
                    return d.alleleIds.map(i => `${slug}_allele_${i}`).join(" ") + ` ${slug}_link`
                  }
                  return `${slug}_link`
                })
                .attr("id", d => getLinkId(d))
                .attr("stroke", "#999999")
                .style("stroke-width", d => `${widthScale(d.percentReads)}em`)
                .style("opacity", d => opacityScale(d.percentReads))

            link.append("title")
              .text(d => {
                if(d.type && d.numReads && d.percentReads) {
                  return `${d.type} in ${d.numReads} (${d.percentReads.toFixed(2)}%) Reads`
                }
                return "Reference"
              })
            return link
          },
        )

  function getAlignmentBounds(vs, c) {
    var os = c.offsets;
    if (c.axis === 'x') {
      var x = vs[os[0].node].x;
      c.bounds = new cola.Rectangle(x, x,
                                    Math.min.apply(Math, os.map(function (o) { return vs[o.node].bounds.y - 20; })),
                                    Math.max.apply(Math, os.map(function (o) { return vs[o.node].bounds.Y + 20; })));
    } else if (c.axis === 'y') {
      var y = vs[os[0].node].y;
      c.bounds = new cola.Rectangle(
        Math.min.apply(Math, os.map(function (o) { return vs[o.node].bounds.x - 20; })),
        Math.max.apply(Math, os.map(function (o) { return vs[o.node].bounds.X + 20; })),
        y, y);
    }
    return c.bounds;
  }

  const createDeletions = (svg, setcolaResult) => svg.selectAll(`.${slug}_deletion`)
        .data(setcolaResult.links.filter(d => d.type === "Deletion"), d => getLinkId(d))
        .join(
          enter => {
            let deletion = enter.append("path")
                .attr("class", d => {
                  if(d.alleleIds) {
                    return d.alleleIds.map(i => `${slug}_allele_${i}`).join(" ") + ` ${slug}_deletion ${slug}_link ${getLinkId(d)}`
                  }
                  return `${slug}_deletion ${slug}_link ${getLinkId(d)}`
                })
                .attr("fill", "none")
                .attr("stroke", "#999999")
                .attr("stroke-dasharray", "5,5")
                .style("stroke-width", d => `${widthScale(d.percentReads)}em`)
                .style("opacity", d => opacityScale(d.percentReads))
            deletion.append("title")
              .text(d => `${d.type} in ${d.numReads} (${d.percentReads.toFixed(2)}%) Reads `)
            return deletion
          }
        )

  const createCleavagePosition = (svg, setcolaResult) => svg.selectAll(`#${slug}_cleavage_position`)
        .data(setcolaResult.nodes.filter(d => d.cleavagePosition))
        .join("line")
        .attr("class", `${slug}_guideline`)
        .attr("stroke-dasharray", "8,8")
        .attr("stroke-width", 1)
        .attr("stroke", "#999999")

  const createCleavageLabel = svg => svg.selectAll(`#${slug}_cleavage_label`)
        .data([{"text": "Predicted cleavage position"}])
        .join("text")
        .attr("id", `${slug}_cleavage_label`)
        .attr("class", `${slug}_group_label`)
        .text(d => d.text)

  const createNodeGroup = (svg, setcolaResult) => svg.selectAll(`.${slug}_node`)
        .data(setcolaResult.nodes, d => d["id"])
        .join(
          enter => {
            let group = enter.append("g")
                .attr("class", d => {
                  if(d.alleleIds) {
                    return d.alleleIds.map(i => `${slug}_allele_${i}`).join(" ") + ` ${slug}_node`
                  }
                  return `${slug}_node ${slug}_reference`
                })
                .attr("id", d => `${slug}_node_${d.id}`)
            group.append("rect")
              .attr("width", d => d.width)
              .attr("height", d => nodeHeight)
              .attr("rx", 5).attr("ry", 5)
              .style("fill", d => color[d.name])
              .style("fill-opacity", d => opacityScale(d.percentReads))
              .call(d3cola.drag)
            group.append("text")
              .attr("class", `${slug}_label`)
              .attr("text-anchor", "middle")
              .style("font-size", "0.5em")
              .text(d => d.name)
              .call(d3cola.drag)
            group.append("title")
              .text(d => {
                if(d.type && d.numReads && d.percentReads) {
                  return `${d.type} Node (${d.id}) in ${d.numReads} (${d.percentReads.toFixed(2)}%) Reads`
                }
                return `Node (${d.id})`
              });
            group
              .style("opacity", d => isInvisibleNode(d) ? 0 : 1)
            return group
          },
          update => {
            update.filter(isInvisibleNode)
              .remove()
            update.select("rect")
              .style("fill", d => d.collapsed ? "#C8C8C8" : color[d.name])
              .style("fill-opacity", d => opacityScale(d.percentReads))
              .attr("width", d => {
                if(d.collapsed && d.collapsedNodes && d.collapsedNodes.length <= 3) {
                  let width = 10 * d.collapsedNodes.length
                  d.width = width
                  return width
                } else {
                  d.width = 15
                  return 15
                }
              })
            update.select("text")
              .text(d => {
                if(d.collapsed) {
                  if(d.collapsedNodes.length > 3) {
                    return d.collapsedNodes.length
                  } else {
                    return d.collapsedNodes.map(e => graph.nodes[e].name).join("")
                  }
                } else {
                  return d.name
                }
              })
            update
              .style("opacity", d => isInvisibleNode(d) ? 0 : 1)
            return update
          },
          exit => exit.remove()
        )

  const nodeHeight = 10

  updateGraph(svg, d3cola)

}

{% for amplicon in report_data['figures']['sgRNA_based_names'] %}
{% for fig_name in report_data['figures']['sgRNA_based_names'][amplicon]['9b'] %}
buildGraphQuilt({{fig_name}}, "{{fig_name}}")
{% endfor %}
{% endfor %}
