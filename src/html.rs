use crate::report::Report;

const UPLOT_CSS: &str = include_str!("../assets/uPlot.min.css");
const UPLOT_JS: &str = include_str!("../assets/uPlot.iife.min.js");

pub fn render_html(report: &Report) -> Result<String, String> {
    let json = serde_json::to_string(report)
        .map_err(|error| format!("failed to serialize report for HTML: {error}"))?
        .replace("</", "<\\/");
    let mut html = String::with_capacity(UPLOT_JS.len() + json.len() + 36_000);
    html.push_str(HTML_HEAD);
    html.push_str(UPLOT_CSS);
    html.push_str("</style></head><body>");
    html.push_str(HTML_BODY);
    html.push_str("<script>");
    html.push_str(UPLOT_JS);
    html.push_str("</script><script>const REPORT=");
    html.push_str(&json);
    html.push(';');
    html.push_str(REPORT_JS);
    html.push_str("</script></body></html>");
    Ok(html)
}

const HTML_HEAD: &str = r#"<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>RustQC report</title>
<style>
:root{color-scheme:light;--ink:#162328;--muted:#637278;--line:#d7e0e2;--paper:#fff;--wash:#edf2f2;--surface:#f8fbfb;--pass:#147a57;--warn:#ae6708;--fail:#c33d3d;--skip:#76848a;--blue:#1f669e;--teal:#087f7a;--violet:#7351a2;--orange:#c05c24}
*{box-sizing:border-box;letter-spacing:0}
html{scroll-behavior:smooth;background:var(--wash)}
body{margin:0;color:var(--ink);font-family:Inter,ui-sans-serif,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;font-size:14px;line-height:1.45}
header{height:92px;background:#12292d;color:#fff;border-bottom:4px solid #d4ae45;display:flex;align-items:center;padding:0 32px;position:relative;z-index:2;box-shadow:0 2px 8px rgba(18,41,45,.12)}
.brand{font-size:25px;font-weight:780}.brand span{color:#e1bf59}.header-file{margin-left:30px;min-width:0;color:#c6d0d1;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}.header-meta{margin-left:auto;color:#c6d0d1;font-variant-numeric:tabular-nums}
.layout{display:grid;grid-template-columns:256px minmax(0,1fr);max-width:1580px;margin:0 auto;background:var(--paper);min-height:calc(100vh - 92px);box-shadow:0 0 0 1px rgba(18,41,45,.04)}
nav{border-right:1px solid var(--line);background:#f7fafa;padding:24px 16px;position:sticky;top:0;height:100vh;overflow:auto}
nav h2{font-size:11px;text-transform:uppercase;color:var(--muted);margin:0 10px 12px;font-weight:760}
nav a{display:grid;grid-template-columns:minmax(0,1fr) auto;gap:10px;align-items:center;padding:8px 10px;color:var(--ink);text-decoration:none;border-left:2px solid transparent;font-size:13px;border-radius:0 4px 4px 0}
nav a:hover{background:#eaf1f1;border-left-color:#547075}.sidebar-status{font-size:9px;line-height:19px;height:19px;min-width:38px;padding:0 6px;box-shadow:none}
main{min-width:0}.overview{padding:34px 42px 30px;background:#fcfefe;border-bottom:1px solid var(--line)}
.overview-kicker{margin:0 0 5px;color:var(--teal);font-size:11px;font-weight:800;text-transform:uppercase}h1{font-size:29px;margin:0 0 5px;font-weight:760}.subtitle{color:var(--muted);margin:0;font-variant-numeric:tabular-nums}
.metrics{display:grid;grid-template-columns:repeat(5,minmax(110px,1fr));margin-top:26px;border:1px solid var(--line);background:#fff;box-shadow:0 3px 12px rgba(18,41,45,.045)}.metric{padding:15px 17px;border-right:1px solid var(--line);min-width:0}.metric:last-child{border-right:0}.metric-label{font-size:11px;color:var(--muted);text-transform:uppercase;font-weight:760}.metric-value{font-size:21px;font-weight:720;margin-top:2px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;font-variant-numeric:tabular-nums}
.module{padding:32px 42px 40px;border-bottom:1px solid var(--line);scroll-margin-top:12px}.module-head{display:flex;align-items:center;gap:12px;margin-bottom:20px}.module h2{font-size:18px;margin:0;font-weight:740}.status{font-size:10px;line-height:22px;height:22px;min-width:48px;padding:0 8px;text-align:center;color:#fff;font-weight:800;border-radius:3px;box-shadow:0 1px 2px rgba(18,41,45,.16)}.status.pass{background:var(--pass)}.status.warn{background:var(--warn)}.status.fail{background:var(--fail)}.status.skip{background:var(--skip)}
.plot{width:100%;height:360px;min-height:360px;position:relative;overflow:hidden;border:1px solid var(--line);background:var(--surface);box-shadow:0 2px 7px rgba(18,41,45,.035)}.plot.tall{height:430px;min-height:430px}.empty{height:120px;display:grid;place-items:center;color:var(--muted);border:1px dashed var(--line);background:var(--surface)}
.chart-legend{display:flex;flex-wrap:wrap;gap:8px 16px;margin:13px 2px 0;color:#44565c;font-size:12px;font-variant-numeric:tabular-nums}.legend-item{display:inline-flex;align-items:center;gap:6px;white-space:nowrap}.legend-swatch{display:block;width:18px;height:3px;border-radius:2px;background:var(--skip);box-shadow:0 0 0 1px rgba(16,41,45,.1)}
.chart-tooltip,.tile-tooltip{position:absolute;z-index:5;min-width:184px;max-width:260px;padding:10px 12px;background:rgba(18,41,45,.96);color:#f7fbfb;box-shadow:0 7px 20px rgba(10,29,31,.22);font-size:12px;line-height:1.35;pointer-events:none}.chart-tooltip[hidden],.tile-tooltip[hidden]{display:none}.tooltip-title{margin-bottom:6px;color:#bfcfce;font-size:11px;font-weight:720}.tooltip-row{display:flex;align-items:center;justify-content:space-between;gap:16px;white-space:nowrap}.tooltip-name{display:flex;align-items:center;gap:6px;min-width:0}.tooltip-marker{width:7px;height:7px;border-radius:50%;flex:0 0 auto}.tooltip-value{font-variant-numeric:tabular-nums;font-weight:700}
table{width:100%;border-collapse:collapse;font-variant-numeric:tabular-nums}th,td{text-align:left;padding:9px 12px;border-bottom:1px solid #e5ebec;vertical-align:top}th{font-size:11px;text-transform:uppercase;color:var(--muted);background:#f4f8f8;position:sticky;top:0}td code{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:12px;overflow-wrap:anywhere}.table-wrap{max-height:430px;overflow:auto;border:1px solid var(--line)}
.basic-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));border-top:1px solid var(--line);border-left:1px solid var(--line)}.basic-row{display:grid;grid-template-columns:minmax(130px,42%) 1fr;border-right:1px solid var(--line);border-bottom:1px solid var(--line)}.basic-key,.basic-value{padding:10px 12px}.basic-key{background:#f4f8f8;color:var(--muted);font-weight:680}.basic-value{overflow-wrap:anywhere}
.tile-wrap{width:100%;position:relative;overflow:auto;border:1px solid var(--line);padding:12px;background:var(--surface)}.tile-wrap canvas{display:block;width:100%;height:420px;background:#fff}.heatmap-legend{display:flex;align-items:center;flex-wrap:wrap;gap:7px;margin:12px 2px 1px;color:#44565c;font-size:12px}.heatmap-label{margin-right:3px;font-weight:700}.heatmap-scale{display:flex;gap:2px}.heatmap-swatch{width:18px;height:9px}.heatmap-low{background:#c33d3d}.heatmap-midlow{background:#e79b8f}.heatmap-mid{background:#e9eeee}.heatmap-midhigh{background:#9ac7d2}.heatmap-high{background:#1f669e}
.uplot{font-family:inherit!important}.u-title,.u-legend{display:none!important}
@media(max-width:900px){header{height:auto;min-height:80px;padding:16px 20px;flex-wrap:wrap}.header-file{order:3;width:100%;margin:6px 0 0}.layout{display:block}.layout>nav{position:relative;height:auto;border-right:0;border-bottom:1px solid var(--line);display:none}.overview,.module{padding-left:20px;padding-right:20px}.metrics{grid-template-columns:repeat(2,minmax(0,1fr))}.metric{border-bottom:1px solid var(--line)}.basic-grid{grid-template-columns:1fr}.plot{height:320px;min-height:320px}.chart-tooltip{max-width:220px}.tile-wrap{padding:8px}}
@media print{header{background:#fff;color:#000;border-bottom:2px solid #000}.header-file,.header-meta{color:#333}.layout{display:block}.layout>nav{display:none}.module{break-inside:avoid}.uplot{max-width:100%}}
"#;

const HTML_BODY: &str = r#"
<header><div class="brand">Rust<span>QC</span></div><div class="header-file" id="header-file"></div><div class="header-meta" id="header-meta"></div></header>
<div class="layout">
<nav><h2>Modules</h2><div id="nav-modules"></div></nav>
<main>
  <section class="overview"><p class="overview-kicker">Analysis summary</p><h1>Sequence quality report</h1><p class="subtitle" id="subtitle"></p><div class="metrics" id="metrics"></div></section>
  <section class="module" id="basic-statistics"><div class="module-head"><span class="status" data-status="basic_statistics"></span><h2>Basic Statistics</h2></div><div class="basic-grid" id="basic-table"></div></section>
  <section class="module" id="per-base-quality"><div class="module-head"><span class="status" data-status="per_base_sequence_quality"></span><h2>Per base sequence quality</h2></div><div class="plot" id="quality-plot"></div></section>
  <section class="module" id="per-tile-quality"><div class="module-head"><span class="status" data-status="per_tile_sequence_quality"></span><h2>Per tile sequence quality</h2></div><div class="tile-wrap" id="tile-wrap"><canvas id="tile-canvas"></canvas></div></section>
  <section class="module" id="per-sequence-quality"><div class="module-head"><span class="status" data-status="per_sequence_quality_scores"></span><h2>Per sequence quality scores</h2></div><div class="plot" id="sequence-quality-plot"></div></section>
  <section class="module" id="per-base-content"><div class="module-head"><span class="status" data-status="per_base_sequence_content"></span><h2>Per base sequence content</h2></div><div class="plot" id="base-content-plot"></div></section>
  <section class="module" id="per-sequence-gc"><div class="module-head"><span class="status" data-status="per_sequence_gc_content"></span><h2>Per sequence GC content</h2></div><div class="plot" id="gc-plot"></div></section>
  <section class="module" id="per-base-n"><div class="module-head"><span class="status" data-status="per_base_n_content"></span><h2>Per base N content</h2></div><div class="plot" id="n-plot"></div></section>
  <section class="module" id="sequence-length"><div class="module-head"><span class="status" data-status="sequence_length_distribution"></span><h2>Sequence Length Distribution</h2></div><div class="plot" id="length-plot"></div></section>
  <section class="module" id="duplication"><div class="module-head"><span class="status" data-status="sequence_duplication_levels"></span><h2>Sequence Duplication Levels</h2></div><div class="plot" id="duplication-plot"></div></section>
  <section class="module" id="overrepresented"><div class="module-head"><span class="status" data-status="overrepresented_sequences"></span><h2>Overrepresented sequences</h2></div><div id="overrepresented-table"></div></section>
  <section class="module" id="adapter-content"><div class="module-head"><span class="status" data-status="adapter_content"></span><h2>Adapter Content</h2></div><div class="plot tall" id="adapter-plot"></div></section>
</main></div>
"#;

const REPORT_JS: &str = r##"
const MODULES=[
  ["basic_statistics","basic-statistics","Basic Statistics"],
  ["per_base_sequence_quality","per-base-quality","Per base sequence quality"],
  ["per_tile_sequence_quality","per-tile-quality","Per tile sequence quality"],
  ["per_sequence_quality_scores","per-sequence-quality","Per sequence quality scores"],
  ["per_base_sequence_content","per-base-content","Per base sequence content"],
  ["per_sequence_gc_content","per-sequence-gc","Per sequence GC content"],
  ["per_base_n_content","per-base-n","Per base N content"],
  ["sequence_length_distribution","sequence-length","Sequence Length Distribution"],
  ["sequence_duplication_levels","duplication","Sequence Duplication Levels"],
  ["overrepresented_sequences","overrepresented","Overrepresented sequences"],
  ["adapter_content","adapter-content","Adapter Content"]
];
const basic=REPORT.basic_statistics.data;
const fmt=new Intl.NumberFormat();
const valueFmt=new Intl.NumberFormat(undefined,{maximumFractionDigits:2});
function localTimestamp(unixSeconds){const date=new Date(unixSeconds*1000);const pad=value=>String(value).padStart(2,"0");const offset=-date.getTimezoneOffset();const sign=offset>=0?"+":"-";const timezone=Intl.DateTimeFormat().resolvedOptions().timeZone||"Local time";return `${date.getFullYear()}-${pad(date.getMonth()+1)}-${pad(date.getDate())}T${pad(date.getHours())}:${pad(date.getMinutes())}:${pad(date.getSeconds())}${sign}${pad(Math.floor(Math.abs(offset)/60))}:${pad(Math.abs(offset)%60)} (${timezone})`}
document.title=`RustQC · ${basic.filename}`;
document.getElementById("header-file").textContent=basic.filename;
document.getElementById("header-meta").textContent=`RustQC ${REPORT.version}`;
document.getElementById("subtitle").textContent=`Generated ${localTimestamp(REPORT.generated_at_unix)}`;

for(const [key,id,label] of MODULES){
  const status=REPORT[key].status;
  const badge=document.querySelector(`[data-status="${key}"]`);
  badge.textContent=status.toUpperCase(); badge.classList.add(status);
  const a=document.createElement("a"); a.href=`#${id}`;
  const text=document.createElement("span"); text.textContent=label;
  const sidebarStatus=document.createElement("span");sidebarStatus.className=`status sidebar-status ${status}`;sidebarStatus.textContent=status.toUpperCase();
  a.append(text,sidebarStatus); document.getElementById("nav-modules").append(a);
}

function metric(label,value){const el=document.createElement("div");el.className="metric";const l=document.createElement("div");l.className="metric-label";l.textContent=label;const v=document.createElement("div");v.className="metric-value";v.textContent=value;el.append(l,v);return el}
const metrics=document.getElementById("metrics");
metrics.append(metric("Sequences",fmt.format(basic.total_sequences)),metric("Total bases",fmt.format(basic.total_bases)),metric("Mean length",basic.mean_length.toFixed(1)),metric("GC",`${basic.gc_percent.toFixed(1)}%`),metric("Encoding",basic.encoding));

function basicRow(key,value){const row=document.createElement("div");row.className="basic-row";const k=document.createElement("div");k.className="basic-key";k.textContent=key;const v=document.createElement("div");v.className="basic-value";v.textContent=value;row.append(k,v);return row}
const basicTable=document.getElementById("basic-table");
for(const [key,value] of [["Filename",basic.filename],["File type",basic.file_type],["Encoding",basic.encoding],["Total Sequences",fmt.format(basic.total_sequences)],["Total Bases",fmt.format(basic.total_bases)],["Sequences flagged as poor quality",fmt.format(basic.sequences_flagged_as_poor_quality)],["Sequence length",basic.sequence_length],["GC content",`${basic.gc_percent.toFixed(2)}%`]]) basicTable.append(basicRow(key,value));

const plots=[];
function formatValue(value){return Number.isFinite(value)?valueFmt.format(value):"No data"}
function appendTooltipRow(host,series,value){const row=document.createElement("div");row.className="tooltip-row";const name=document.createElement("span");name.className="tooltip-name";const marker=document.createElement("span");marker.className="tooltip-marker";marker.style.background=series.color;const label=document.createElement("span");label.textContent=series.label;name.append(marker,label);const text=document.createElement("span");text.className="tooltip-value";text.textContent=formatValue(value);row.append(name,text);host.append(row)}
function updateChartTooltip(plot,el,tooltip,labels,series){const index=plot.cursor.idx;if(index==null){tooltip.hidden=true;return}const title=document.createElement("div");title.className="tooltip-title";title.textContent=labels[index]??`Position ${formatValue(plot.data[0][index])}`;tooltip.replaceChildren(title);series.forEach((series,index)=>appendTooltipRow(tooltip,series,plot.data[index+1][plot.cursor.idx]));const width=tooltip.offsetWidth||220;const height=tooltip.offsetHeight||120;const left=Math.max(10,Math.min(el.clientWidth-width-10,(plot.cursor.left||0)+18));const top=Math.max(10,Math.min(el.clientHeight-height-10,(plot.cursor.top||0)+18));tooltip.style.left=`${left}px`;tooltip.style.top=`${top}px`;tooltip.hidden=false}
function chartLegend(el,series){const legend=document.createElement("div");legend.className="chart-legend";legend.setAttribute("aria-label","Chart legend");for(const seriesItem of series){const item=document.createElement("span");item.className="legend-item";const swatch=document.createElement("span");swatch.className="legend-swatch";swatch.style.background=seriesItem.color;const label=document.createElement("span");label.textContent=seriesItem.label;item.append(swatch,label);legend.append(item)}el.insertAdjacentElement("afterend",legend)}
function linePlot(id,xLabel,x,series,yLabel,range,labels){
  const el=document.getElementById(id);if(!x.length){el.className="empty";el.textContent="No observations";return}
  const tooltip=document.createElement("div");tooltip.className="chart-tooltip";tooltip.hidden=true;tooltip.setAttribute("role","status");el.append(tooltip);
  const options={width:Math.max(320,el.clientWidth),height:el.clientHeight,scales:{x:{time:false},y:range?{range}:{}},axes:[{label:xLabel},{label:yLabel}],series:[{},...series.map(s=>({label:s.label,stroke:s.color,width:s.width||2,dash:s.dash||[],points:{show:false}}))],legend:{show:false},cursor:{drag:{x:true,y:false}},hooks:{setCursor:[plot=>updateChartTooltip(plot,el,tooltip,labels,series)]}};
  const plot=new uPlot(options,[x,...series.map(s=>s.values)],el);plots.push([plot,el]);chartLegend(el,series);
}
function histogram(id,points,xLabel,yLabel,color){linePlot(id,xLabel,points.map(p=>p.value),[{label:yLabel,values:points.map(p=>p.count),color}],yLabel,undefined,points.map(p=>p.label))}

const quality=REPORT.per_base_sequence_quality.data;
linePlot("quality-plot","Base position",quality.map(p=>p.position),[
  {label:"Mean",values:quality.map(p=>p.mean),color:"#2469a0",width:2.5},
  {label:"Median",values:quality.map(p=>p.median),color:"#16835b"},
  {label:"Lower quartile",values:quality.map(p=>p.lower_quartile),color:"#c65f25"},
  {label:"Upper quartile",values:quality.map(p=>p.upper_quartile),color:"#7a56a6"},
  {label:"10th percentile",values:quality.map(p=>p.percentile_10),color:"#c83a3a",dash:[6,4]},
  {label:"90th percentile",values:quality.map(p=>p.percentile_90),color:"#16817b",dash:[6,4]}
],"Phred quality",()=>[0,Math.max(42,...quality.map(p=>p.percentile_90||0))],quality.map(p=>p.base));

const seqQuality=REPORT.per_sequence_quality_scores.data.distribution;
histogram("sequence-quality-plot",seqQuality,"Quality score","Reads","#7a56a6");

const content=REPORT.per_base_sequence_content.data;
linePlot("base-content-plot","Base position",content.map(p=>p.position),[
  {label:"G",values:content.map(p=>p.g),color:"#16835b"},{label:"A",values:content.map(p=>p.a),color:"#2469a0"},{label:"T",values:content.map(p=>p.t),color:"#c83a3a"},{label:"C",values:content.map(p=>p.c),color:"#c65f25"}
],"Percent",()=>[0,100],content.map(p=>p.base));

const gc=REPORT.per_sequence_gc_content.data.distribution;
linePlot("gc-plot","GC content (%)",gc.map(p=>p.gc_percent),[{label:"Observed",values:gc.map(p=>p.count),color:"#2469a0",width:2.5},{label:"Theoretical",values:gc.map(p=>p.theoretical),color:"#c65f25",dash:[7,4]}],"Reads",undefined,gc.map(p=>`${p.gc_percent}%`));

const n=REPORT.per_base_n_content.data;
linePlot("n-plot","Base position",n.map(p=>p.position),[{label:"N",values:n.map(p=>p.value),color:"#7a56a6",width:2.5}],"Percent",()=>[0,Math.max(5,...n.map(p=>p.value))*1.1],n.map(p=>p.base));

histogram("length-plot",REPORT.sequence_length_distribution.data,"Sequence length","Reads","#16817b");

const dup=REPORT.sequence_duplication_levels.data.levels;
linePlot("duplication-plot","Duplication band",dup.map((_,i)=>i+1),[{label:"Percent of total",values:dup.map(p=>p.percentage_of_total),color:"#b66b08",width:2.5}],"Percent",undefined,dup.map(p=>p.level));

const adapters=REPORT.adapter_content.data;
linePlot("adapter-plot","Base position",adapters.position_midpoints,adapters.series.map((s,i)=>({label:s.name,values:s.percentages,color:["#2469a0","#c83a3a","#16835b","#7a56a6","#c65f25","#16817b"][i%6]})),"Percent",()=>[0,Math.max(5,adapters.max_content*1.15)],adapters.position_midpoints.map(value=>`Base ${value}`));

function renderOverrepresented(){const rows=REPORT.overrepresented_sequences.data;const host=document.getElementById("overrepresented-table");if(!rows.length){host.className="empty";host.textContent="No overrepresented sequences";return}const wrap=document.createElement("div");wrap.className="table-wrap";const table=document.createElement("table");const head=document.createElement("thead");const hr=document.createElement("tr");for(const label of ["Sequence","Count","Percentage","Possible source"]){const th=document.createElement("th");th.textContent=label;hr.append(th)}head.append(hr);table.append(head);const body=document.createElement("tbody");for(const row of rows){const tr=document.createElement("tr");const seq=document.createElement("td");const code=document.createElement("code");code.textContent=row.sequence;seq.append(code);const count=document.createElement("td");count.textContent=fmt.format(row.count);const pct=document.createElement("td");pct.textContent=`${row.percentage.toFixed(3)}%`;const source=document.createElement("td");source.textContent=row.possible_source;tr.append(seq,count,pct,source);body.append(tr)}table.append(body);wrap.append(table);host.append(wrap)}
renderOverrepresented();

function heatmapLegend(wrap){const old=wrap.querySelector(".heatmap-legend");if(old)old.remove();const legend=document.createElement("div");legend.className="heatmap-legend";const label=document.createElement("span");label.className="heatmap-label";label.textContent="Mean quality deviation";const low=document.createElement("span");low.textContent="Lower";const scale=document.createElement("span");scale.className="heatmap-scale";for(const name of ["heatmap-low","heatmap-midlow","heatmap-mid","heatmap-midhigh","heatmap-high"]){const swatch=document.createElement("span");swatch.className=`heatmap-swatch ${name}`;scale.append(swatch)}const high=document.createElement("span");high.textContent="Higher";legend.append(label,low,scale,high);wrap.append(legend)}
function updateTileTooltip(event){const state=window.rustqcTileState;if(!state)return;const {canvas,tooltip,module,left,top,cellW,cellH,width,height}=state;const rect=canvas.getBoundingClientRect();const x=(event.clientX-rect.left)*width/rect.width;const y=(event.clientY-rect.top)*height/rect.height;const column=Math.floor((x-left)/cellW);const row=Math.floor((y-top)/cellH);if(column<0||column>=module.data.positions.length||row<0||row>=module.data.tiles.length){tooltip.hidden=true;return}const value=module.data.deviations[row][column];const title=document.createElement("div");title.className="tooltip-title";title.textContent=`Tile ${module.data.tiles[row]} | Base ${module.data.positions[column]}`;const detail=document.createElement("div");detail.textContent=`Deviation: ${value.toFixed(2)} Phred`;tooltip.replaceChildren(title,detail);const tooltipWidth=tooltip.offsetWidth||210;const tooltipHeight=tooltip.offsetHeight||56;const tooltipLeft=Math.max(8,Math.min(rect.width-tooltipWidth-8,event.clientX-rect.left+14));const tooltipTop=Math.max(8,Math.min(rect.height-tooltipHeight-8,event.clientY-rect.top+14));tooltip.style.left=`${tooltipLeft}px`;tooltip.style.top=`${tooltipTop}px`;tooltip.hidden=false}
function renderTile(){const module=REPORT.per_tile_sequence_quality;const wrap=document.getElementById("tile-wrap");if(module.status==="skip"||!module.data.tiles.length){wrap.className="empty";wrap.textContent="Tile identifiers unavailable";return}wrap.className="tile-wrap";const canvas=document.getElementById("tile-canvas");const rect=canvas.getBoundingClientRect();const dpr=window.devicePixelRatio||1;const width=Math.max(720,rect.width);const height=420;canvas.width=width*dpr;canvas.height=height*dpr;const ctx=canvas.getContext("2d");ctx.scale(dpr,dpr);const left=56,top=12,right=10,bottom=28;const cellW=(width-left-right)/module.data.positions.length;const cellH=(height-top-bottom)/module.data.tiles.length;const max=Math.max(5,module.data.max_deviation);for(let y=0;y<module.data.tiles.length;y++){for(let x=0;x<module.data.positions.length;x++){const v=module.data.deviations[y][x];const strength=Math.min(1,Math.abs(v)/max);ctx.fillStyle=v<0?`rgba(195,61,61,${0.12+strength*.88})`:`rgba(31,102,158,${0.12+strength*.88})`;ctx.fillRect(left+x*cellW,top+y*cellH,Math.ceil(cellW),Math.ceil(cellH))}}ctx.fillStyle="#4d5e63";ctx.font="11px sans-serif";ctx.textAlign="right";const tick=Math.max(1,Math.ceil(module.data.tiles.length/12));for(let y=0;y<module.data.tiles.length;y+=tick)ctx.fillText(module.data.tiles[y],left-6,top+(y+.75)*cellH);ctx.textAlign="center";const xTick=Math.max(1,Math.ceil(module.data.positions.length/10));for(let x=0;x<module.data.positions.length;x+=xTick)ctx.fillText(module.data.positions[x],left+(x+.5)*cellW,height-8);heatmapLegend(wrap);let tooltip=wrap.querySelector(".tile-tooltip");if(!tooltip){tooltip=document.createElement("div");tooltip.className="tile-tooltip";tooltip.hidden=true;tooltip.setAttribute("role","status");wrap.append(tooltip)}window.rustqcTileState={canvas,tooltip,module,left,top,cellW,cellH,width,height};canvas.onmousemove=updateTileTooltip;canvas.onmouseleave=()=>{tooltip.hidden=true}}
renderTile();

let resizeTimer;addEventListener("resize",()=>{clearTimeout(resizeTimer);resizeTimer=setTimeout(()=>{for(const [plot,el] of plots)plot.setSize({width:Math.max(320,el.clientWidth),height:el.clientHeight});renderTile()},120)});
"##;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn embedded_assets_are_present() {
        assert!(UPLOT_JS.len() > 10_000);
        assert!(UPLOT_CSS.contains(".uplot"));
        assert!(REPORT_JS.contains("chart-legend"));
        assert!(REPORT_JS.contains("chart-tooltip"));
        assert!(REPORT_JS.contains("heatmap-legend"));
        assert!(REPORT_JS.contains("localTimestamp"));
        assert!(REPORT_JS.contains("sidebar-status"));
    }
}
