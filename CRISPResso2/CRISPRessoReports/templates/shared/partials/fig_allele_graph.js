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

  const seqLength = graph["alleles"][0].seq.length,
        tableCellHeight = 10,
        tableCellWidth = 10,
        numGuides = graph["groups"].length
  const guidePadding = tableCellHeight / 2,
        sliderWidth = tableCellWidth * 2
  const guideHeight = numGuides > 0 ? (guidePadding * 2) + (numGuides * tableCellHeight) : 0
  const quiltHeight = (graph["alleles"].length * tableCellHeight) + guideHeight

  const quiltSvg = d3.select(tableId)
        .append("svg")
        .attr("viewBox", `0 -10 ${(seqLength * tableCellWidth) + sliderWidth + 170} ${quiltHeight}`)
        .style("user-select", "none")

  const alleleGroup = quiltSvg.selectAll(`.${slug}_allele`)
        .data(graph["alleles"], d => d.id)
        .join(
          enter => {
            let root = enter.append("g")
            root.style("opacity", 1)
              .attr("transform", (d, i) => {
                if(i == 0) {
                  return `translate(${sliderWidth},0)`
                } else {
                  return `translate(${sliderWidth},${(i * tableCellHeight) + guideHeight})`
                }
              })
              .attr("id", d => `${slug}_allele_${d.id}`)
            root.append("text")
              .attr("x", (seqLength * tableCellWidth) + 2)
              .attr("dy", "-0.2em")
              .style("font-family", "Arial")
              .style("font-size", "7pt")
              .style("fill", "#262626")
              .text(d => `${d.name}, ${d.numReads} (${d.percentReads.toFixed(2)}%) Reads`)
            return root
          },
        )
  const guides = quiltSvg.selectAll(`.${slug}_guide`)
        .data(graph["sgRNAs"], d => d)
        .join(
          enter => {
            let root = enter.append("g")
            root
              .attr("transform", (d, i) => `translate(${(d["interval"][0] * tableCellWidth) + sliderWidth},${(d["row"] * tableCellHeight) + guidePadding})`)
            root.append("rect")
              .attr("width", d => (d["interval"][1] - d["interval"][0]) * tableCellWidth)
              .attr("height", tableCellHeight)
              .style("fill", "#D9D9D9")
            root.append("text")
              .attr("dx", "0.2em")
              .attr("dy", "1em")
              .attr("text-anchor", "start")
              .style("font-size", "6pt")
              .text(d => d.hasOwnProperty("name") ? d["name"] : "sgRNA")
            return root
          }
        )
        .selectAll(`.${slug}_guide_mismatch`)
        .data(d => d["mismatch"], d => d)
        .join(
          enter => enter.append("rect")
            .attr("transform", d => `translate(${d * tableCellWidth},0)`)
            .attr("width", tableCellWidth)
            .attr("height", tableCellHeight)
            .style("fill", color["highlight"])
            .style("opacity", "0.5")
        )
  const alleleSeqGroup = alleleGroup.selectAll(`.${slug}_alleleSeq`)
        .data(d => d.seq, d => d)
        .join(enter => {
          let root = enter.append("g")
          root.append("rect")
            .attr("fill", d => color[d])
            .attr("width", tableCellWidth)
            .attr("height", tableCellHeight)
            .attr("x", (d, i) => i * tableCellWidth)
            .attr("y", -10)
            .style("opacity", "0.5")
          root.append("text")
            .attr("x", (d, i) => (i * tableCellWidth) + (tableCellWidth / 2))
            .attr("text-anchor", "middle")
            .attr("dy", "-0.2em")
            .style("font-size", "6pt")
            .style("font-family", "Arial")
            .style("stroke-width", "0.3pt")
            .style("fill", "#262626")
            .text(d => d)
          return root
        })
  const quiltCutPoints = quiltSvg.selectAll(`.${slug}_quilt_cut_points`)
        .data(graph["cut_points"], d => d)
        .join(
          enter => {
            let group = enter.append("g")
                .attr("class", `${slug}_cut_point`)
            group.append("line")
              .attr("stroke-dasharray", "4,4")
              .attr("stroke-width", 1)
              .attr("stroke", "#999999")
              .attr("x1", d => (d.position * tableCellWidth) + sliderWidth)
              .attr("x2", d => (d.position * tableCellWidth) + sliderWidth)
              .attr("y1", -10)
              .attr("y2", quiltHeight)

            group.append("title")
              .text(d => d.name)
            /* .attr("dx", d => ((d.position) * tableCellWidth) + sliderWidth + (tableCellWidth / 2))
             * .attr("y", guidePadding + tableCellHeight)
             * .style("font-size", "0.5em") */
            return group
          }
        )

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
                              .attr("width", e[1] * tableCellWidth)
                              .attr("height", tableCellHeight)
                              .attr("x", ((e[0] + 1) * tableCellWidth) + sliderWidth)
                              .attr("y", ((d.id - 1) * tableCellHeight) + guideHeight)
                              .attr("class", `${slug}_allele_${d.id}_insertion`)
                              .style("stroke", color["highlight"])
                              .style("fill", "none")
                             )
    })

  let inactiveAlleles = new Set()

  const updateAlleleGroup = (alleleId, alleleElements, newOpacity, show) => {
    d3.select(`#${slug}_allele_${alleleId}`).style("opacity", newOpacity)
    d3.select(`.${slug}_allele_${alleleId}_insertion`).style("opacity", newOpacity)
    alleleElements.forEach(d => {
      if(d.hasOwnProperty("type") && d["type"] === "Deletion") {
        graph.links[d._linkid].show = show
      }
      else if(d.hasOwnProperty("id")) {
        graph.nodes[d.id].show = show
      }
    })
  }
  const hideAlleleGroup = alleleId => {
    let alleleElements = d3.select(graphId)
        .selectAll(`.${slug}_allele_${alleleId}`)
        .data()
        .filter(d =>
          d.alleleIds.every(e => inactiveAlleles.has(e))
        )

    updateAlleleGroup(alleleId, alleleElements, 0.3, false)
  }
  const showAlleleGroup = alleleId => {
    let alleleElements = d3.select(graphId)
        .selectAll(`.${slug}_allele_${alleleId}`)
        .data()
        .filter(d => {
          if(d.hasOwnProperty("type") && d["type"] === "Deletion") {
            return !graph.links[d._linkid].show
          }
          else if(d.hasOwnProperty("id")) {
            return !graph.nodes[d.id].show
          }
          return false
        })
    updateAlleleGroup(alleleId, alleleElements, 1, true)
  }
  const hideOrShowAlleleGroup = function(d, el, update = false) {
    // check if this allele is currently hidden
    if(!inactiveAlleles.has(d.id)) {
      inactiveAlleles.add(d.id)
      hideAlleleGroup(d.id)
    }
    else {
      inactiveAlleles.delete(d.id)
      // filter the elements if they are already hidden
      showAlleleGroup(d.id)
    }
    if(update) {
      updateGraph(svg, d3cola)
    }
  }

  const setupSlider = (v1, v2) => {
    var sliderVals = [v1, v2]
    const height = (graph["alleles"].length - 2) * tableCellHeight

    var x = d3.scaleLinear()
        .domain([1, graph["alleles"].length])
        .range([0, height])
        .clamp(true);

    var xMin=x(1),
        xMax=x(graph["alleles"].length)

    var slider = quiltSvg.append("g")
        .attr("class", "slider")
        .attr("transform", `translate(${tableCellWidth / 2},${guideHeight}) rotate(90)`);

    slider.append("line")
      .attr("class", "track")
      .attr("x1", x.range()[0] + 5)
      .attr("x2", x.range()[1] + 5)
      .style("stroke", "#DDDDDD")
      .style("stroke-width", 8)
      .style("stroke-linecap", "round")

    var selRange = slider.append("line")
        .attr("class", "sel-range")
        .attr("x1", 10 + x(sliderVals[0]))
        .attr("x2", x(sliderVals[1]))
        .style("stroke", "#0D6EFD")
        .style("stroke-width", 8)

    var handle = slider.selectAll("rect")
        .data([0, 1])
        .enter().append("rect", ".track-overlay")
        .attr("class", "handle")
        .attr("y", -5)
        .attr("x", d => x(sliderVals[d]))
        .attr("rx", 3)
        .attr("height", 10)
        .attr("width", 10)
        .style("fill", "#BBBBBB")
        .call(
          d3.drag()
            .on("start", startDrag)
            .on("drag", drag)
            .on("end", endDrag)
        );

    function startDrag(){
      d3.select(this).raise().classed("active", true);
    }

    function drag(d){
      var x1 = d3.event.x
      if(x1 > xMax){
        x1 = xMax
      }
      else if(x1 < xMin){
        x1 = xMin
      }
      d3.select(this).attr("x", x1)
      var x2 = x(sliderVals[d == 0 ? 1 : 0])
      selRange
        .attr("x1", 10 + x1)
        .attr("x2", 10 + x2)
    }

    function endDrag(d){
      let v = Math.round(x.invert(d3.event.x)),
          elem = d3.select(this)
      sliderVals[d] = v
      let v1 = Math.min(sliderVals[0], sliderVals[1]),
          v2 = Math.max(sliderVals[0], sliderVals[1])
      elem.classed("active", false)
        .attr("x", Math.ceil(x(v) / 10) * 10);
      selRange
        .attr("x1", 10 + x(v1))
        .attr("x2", 10 + x(v2))

      updateAllelesFromSlider(v1, v2)
      updateGraph(svg, d3cola)
    }

    const updateAllelesFromSlider = (start, end) => {
      for(let alleleId = x.domain()[0]; alleleId < x.domain()[1]; alleleId++) {
        let alleleInRange = alleleId >= start && alleleId <= end
        if(!alleleInRange && !inactiveAlleles.has(alleleId)) {
          inactiveAlleles.add(alleleId)
          hideAlleleGroup(alleleId)
        }
        else if(alleleInRange && inactiveAlleles.has(alleleId)) {
          inactiveAlleles.delete(alleleId)
          showAlleleGroup(alleleId)
        }
      }
    }
    updateAllelesFromSlider(Math.min(sliderVals[0], sliderVals[1]), Math.max(sliderVals[0], sliderVals[1]))
  }

  const thisGraph = d3.select(graphId)
  const alleleInsertions = Object.fromEntries(graph["alleles"].slice(1).map(d => [d.id, new Set(d.insertions.map(e => e[0]))]))
  alleleInsertions[0] = new Set()
  const hoverAlleleGroup = (d, nodeStrokeColor, linkStrokeColor) => {
    if(!inactiveAlleles.has(d.id)) {
      thisGraph
        .selectAll(d.alleleNodeIds.map(e => `#${slug}_node_${e}`).join(","))
        .select("rect")
        .style("stroke", nodeStrokeColor)
      thisGraph
        .selectAll(`.${slug}_link.${slug}_allele_${d.id}`)
        .style("stroke", linkStrokeColor)
      d3.pairs(d.alleleNodeIds).forEach(([sourceNodeId, targetNodeId]) => {
        if(!alleleInsertions[d.id].has(sourceNodeId)) {
          thisGraph
            .select(`#${slug}_link_${sourceNodeId}_${targetNodeId}_Reference`)
            .style("stroke", linkStrokeColor)
        }
      })
    }
  }

  alleleGroup
    .filter(d => d.name != "Reference")
    .on("click", function(d) { hideOrShowAlleleGroup(d, this, true) })
  alleleGroup
    .on("mouseover", d => { hoverAlleleGroup(d, color["highlight"], color["highlight"]) })
    .on("mouseout", d => { hoverAlleleGroup(d, "#ffffff", "#999999") })

  const hideOrShowAllAlleleGroups = () => {
    alleleGroup
      .filter(d => d.name != "Reference")
      .each(function(d) { hideOrShowAlleleGroup(d, this, false) })

    updateGraph(svg, d3cola)
  }

  const width = seqLength * 40,
        height = seqLength * 5,
        nodeHeight = 10,
        groupPadding = nodeHeight * 2

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
      return link.show && isVisibleNode(graph.nodes[link.source]) && isVisibleNode(graph.nodes[link.target])
    } else {
      return link.show && isVisibleNode(graph.nodes[link.source.id]) && isVisibleNode(graph.nodes[link.target.id])
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
          .links(graph.links)
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

  const bounds = (value) => Math.min(Math.max(60, value), Math.max(width, height))

  const updateGraph = (svg, d3cola) => {
    let setcolaResult = updateSetcola(d3cola)

    let group = createGroup(svg, d3cola, setcolaResult)
    let groupLabel = createGroupLabel(svg, setcolaResult)
    let link = createLink(svg, setcolaResult)

    let deletions = createDeletions(svg, setcolaResult)

    let cutPoints = createCutPoints(svg)

    let nodeGroup = createNodeGroup(svg, setcolaResult)
    nodeGroup.raise()

    nodeGroup.filter(d => d.collapsedNodes)
      .on("click", d => {collapseNode(d, true)})

    d3cola.on("tick", () => {
      // this will ensure that the reference nodes stay in the middle of the viewport
      graph.nodes = graph.nodes.map(d => {
        if(d.type === "Reference") {
          d.y = height / 2
        }
        return d
      })

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

      group
        .attr("x", d => {
          thisGraph
            .select(`#${slug}_${d.name}_label`)
            .attr("x", bounds(d.bounds.x))
          return bounds(d.bounds.x) - groupPadding / 4
        })
        .attr("width", d => d.bounds.width() + groupPadding / 2)
        .attr("height", d => d.bounds.height() + groupPadding / 2)

      nodeGroup.select("text")
        .attr("x", d => bounds(d.x))
        .attr("y", d => bounds(d.y) + nodeHeight / 4)

      if(cutPoints) {
        cutPoints
          .attr("transform", d => `translate(${setcolaResult.nodes[d.position].x},0)`)
      }
    });
  }

  const getLinkId = link => {
    if(typeof(link.source) == "number" && typeof(link.target) == "number") {
      return `${slug}_link_${link.source}_${link.target}_${link.type}`
    } else {
      return `${slug}_link_${link.source.id}_${link.target.id}_${link.type}`
    }
  }

  const createGroup = (svg, d3cola, setcolaResult) => svg.selectAll(`.${slug}_group`)
        .data(setcolaResult.groups)
        .join("rect")
        .attr("rx", 8)
        .attr("ry", 8)
        .attr("y", (height / 2) - (groupPadding / 4))
        .attr("class", `${slug}_group`)
        .attr("cursor", "move")
        .style("fill", "#888888")
        .style("opacity", 0.5)
        .call(d3cola.drag)

  const createGroupLabel = (svg, setcolaResult) => svg.selectAll(`.${slug}_groupLabel`)
        .data(setcolaResult.groups)
        .join("text")
        .attr("id", d => `${slug}_${d.name}_label`)
        .attr("class", `${slug}_groupLabel`)
        .attr("y", (d, i) => (height / 2) + groupPadding + 5 + (15 * i))
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
                .style("opacity", d => isInvisibleLink(d) ? 0 : opacityScale(d.percentReads))

            link.append("title")
              .text(d => {
                if(d.type && d.numReads && d.percentReads) {
                  return `${d.type} in ${d.numReads} (${d.percentReads.toFixed(2)}%) Reads`
                }
                return "Reference"
              })
            return link
          },
          update => update
            .style("opacity", d => isInvisibleLink(d) ? 0 : opacityScale(d.percentReads))
        )

  const createDeletions = (svg, setcolaResult) => svg.selectAll(`.${slug}_deletion`)
        .data(setcolaResult.links.filter(d => d.type === "Deletion"), d => getLinkId(d))
        .join(
          enter => {
            let deletion = enter.append("path")
                .attr("id", d => getLinkId(d))
                .attr("class", d => {
                  if(d.alleleIds) {
                    return d.alleleIds.map(i => `${slug}_allele_${i}`).join(" ") + ` ${slug}_deletion ${slug}_link`
                  }
                  return `${slug}_deletion ${slug}_link`
                })
                .attr("fill", "none")
                .attr("stroke", "#999999")
                .attr("stroke-dasharray", "5,5")
                .style("stroke-width", d => `${widthScale(d.percentReads)}em`)
                .style("opacity", d => isInvisibleLink(d) ? 0 : opacityScale(d.percentReads))
            deletion.append("title")
              .text(d => `${d.type} in ${d.numReads} (${d.percentReads.toFixed(2)}%) Reads `)
            return deletion
          },
          update => update
            .style("opacity", d => isInvisibleLink(d) ? 0 : opacityScale(d.percentReads)),
        )

  const createCutPoints = svg => svg.selectAll(`.${slug}_cut_point`)
        .data(graph["cut_points"], d => d)
        .join(
          enter => {
            let group = enter.append("g")
                .attr("class", `${slug}_cut_point`)
            group.append("line")
              .attr("stroke-dasharray", "8,8")
              .attr("stroke-width", 1)
              .attr("stroke", "#999999")
              .attr("x1", -15)
              .attr("x2", -15)
              .attr("y1", -10)
              .attr("y2", height)
            group.append("text")
              .attr("dx", 20)
              .attr("y", height - 10)
              .text(d => d.name)
            return group
          }
        )

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
                .attr("cursor", "move")
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

  updateGraph(svg, d3cola)
  setupSlider(0, 3)
  updateGraph(svg, d3cola)

}

{% for amplicon in report_data['figures']['sgRNA_based_names'] %}
{% for fig_name in report_data['figures']['sgRNA_based_names'][amplicon]['9b'] %}
buildGraphQuilt({{fig_name}}, "{{fig_name}}")
{% endfor %}
{% endfor %}
