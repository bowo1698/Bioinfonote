local function html(text)
  return pandoc.RawBlock('html', text)
end

local function escape_html(text)
  return text:gsub('&', '&amp;'):gsub('<', '&lt;'):gsub('>', '&gt;'):gsub('"', '&quot;'):gsub("'", '&#39;')
end

function CodeBlock(el)
  if el.classes:includes('text') and not el.classes:includes('terminal-output') then
    el.classes:insert('terminal-output')
  end
  return el
end

function Pandoc(doc)
  local input = quarto.doc.input_file
  if not input:match('^docs/') and not input:match('/docs/') then
    return doc
  end

  local title = escape_html(pandoc.utils.stringify(doc.meta.title))
  local top = html('<style>#title-block-header { display: none; }</style><header class="article-frame"><nav class="article-frame-nav" aria-label="Navigasi artikel"><a href="../index.html">Beranda</a><a href="../materi.html">Materi</a></nav><div class="article-frame-title"><h1>' .. title .. '</h1></div><p class="article-frame-byline">Oleh <a href="bio.html">Agus Wibowo</a></p></header>')
  local bottom = html('<nav class="article-frame-footer" aria-label="Navigasi artikel"><a href="../index.html">Beranda</a><a href="../materi.html">Materi</a></nav>')

  table.insert(doc.blocks, 1, top)
  table.insert(doc.blocks, bottom)
  return doc
end
