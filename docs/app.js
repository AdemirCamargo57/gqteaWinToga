/* gqteaWinToga site: latest-release download, tabs, and the hero molecule. */

/* ---------- Latest release: point the download button at the installer zip ---------- */
(function () {
  var REPO = "AdemirCamargo57/gqteaWinToga";
  var btn = document.getElementById("download-btn");
  var label = document.getElementById("download-label");
  var meta = document.getElementById("download-meta");
  if (!btn || !window.fetch) return;

  function sizeMB(bytes) { return (bytes / 1048576).toFixed(0) + " MB"; }
  function niceDate(iso) {
    try {
      return new Date(iso).toLocaleDateString("en-GB", { day: "numeric", month: "long", year: "numeric" });
    } catch (e) { return ""; }
  }

  fetch("https://api.github.com/repos/" + REPO + "/releases/latest", {
    headers: { Accept: "application/vnd.github+json" }
  })
    .then(function (r) { if (!r.ok) throw new Error(r.status); return r.json(); })
    .then(function (rel) {
      var assets = rel.assets || [];
      var asset = assets.find(function (a) { return /INSTALLER\.zip$/i.test(a.name); }) ||
                  assets.find(function (a) { return /\.zip$/i.test(a.name); });
      if (!asset) return;
      var version = (rel.tag_name || "").replace(/^v/i, "");
      btn.href = asset.browser_download_url;
      label.textContent = "Download gqteaWinToga " + version;
      meta.textContent = "Windows 10 or 11, 64-bit. " + sizeMB(asset.size) +
        " zip, released " + niceDate(rel.published_at) + ". Free and open source, MIT license.";
    })
    .catch(function () { /* keep the fallback link to the releases page */ });
})();

/* ---------- Tabs (WAI-ARIA tabs pattern) ---------- */
(function () {
  var root = document.querySelector("[data-tabs]");
  if (!root) return;
  var tabs = Array.prototype.slice.call(root.querySelectorAll('[role="tab"]'));

  function select(tab, focus) {
    tabs.forEach(function (t) {
      var on = t === tab;
      t.setAttribute("aria-selected", on ? "true" : "false");
      t.tabIndex = on ? 0 : -1;
      document.getElementById(t.getAttribute("aria-controls")).hidden = !on;
    });
    if (focus) tab.focus();
  }

  tabs.forEach(function (tab, i) {
    tab.addEventListener("click", function () { select(tab, false); });
    tab.addEventListener("keydown", function (e) {
      var j = null;
      if (e.key === "ArrowRight") j = (i + 1) % tabs.length;
      else if (e.key === "ArrowLeft") j = (i - 1 + tabs.length) % tabs.length;
      else if (e.key === "Home") j = 0;
      else if (e.key === "End") j = tabs.length - 1;
      if (j !== null) { e.preventDefault(); select(tabs[j], true); }
    });
  });
})();

/* ---------- Hero molecule: glycine + two H-bonded waters, ball and stick ---------- */
(function () {
  var canvas = document.getElementById("molecule");
  if (!canvas || !canvas.getContext) return;
  var ctx = canvas.getContext("2d");

  // Element, x, y, z in Angstrom (an XYZ frame, as the program reads them)
  var atoms = [
    ["C",  0.00,  0.00,  0.00],
    ["C",  1.52,  0.00,  0.00],
    ["O",  2.17,  1.04,  0.00],
    ["O",  2.10, -1.21,  0.00],
    ["H",  3.06, -1.10,  0.00],
    ["N", -0.48, -1.38,  0.00],
    ["H", -1.49, -1.37,  0.10],
    ["H", -0.17, -1.85,  0.84],
    ["H", -0.37,  0.52,  0.89],
    ["H", -0.37,  0.52, -0.89],
    // water 1, donating to the carbonyl O
    ["O",  2.75,  3.85,  0.35],
    ["H",  2.55,  2.90,  0.25],
    ["H",  3.68,  4.05,  0.60],
    // water 2, accepting from an N-H
    ["O", -3.30, -1.45,  0.30],
    ["H", -3.75, -0.62,  0.55],
    ["H", -3.80, -2.15,  0.82]
  ];
  var hbonds = [[11, 2], [6, 13]];
  var style = {
    C: { r: 0.42, c: "#5d6873" },
    H: { r: 0.26, c: "#f4f6f7" },
    O: { r: 0.40, c: "#d23b2e" },
    N: { r: 0.41, c: "#3157c9" }
  };
  var cov = { C: 0.76, H: 0.31, O: 0.66, N: 0.71 };

  // center the coordinates
  function mid(k) {
    var v = atoms.map(function (a) { return a[k]; });
    return (Math.min.apply(null, v) + Math.max.apply(null, v)) / 2;
  }
  var cx = mid(1), cy = mid(2), cz = mid(3);
  var P = atoms.map(function (a) { return { el: a[0], x: a[1] - cx, y: a[2] - cy, z: a[3] - cz }; });

  // fit the model in the canvas at any orientation
  var maxR = 0;
  P.forEach(function (p) { maxR = Math.max(maxR, Math.sqrt(p.x * p.x + p.y * p.y + p.z * p.z) + 0.45); });

  var CAM = 30; // camera distance for a gentle perspective

  // covalent bonds by distance, like a connection cutoff
  var bonds = [];
  for (var i = 0; i < P.length; i++) for (var j = i + 1; j < P.length; j++) {
    var dx = P[i].x - P[j].x, dy = P[i].y - P[j].y, dz = P[i].z - P[j].z;
    var d = Math.sqrt(dx * dx + dy * dy + dz * dz);
    if (d < 1.2 * (cov[P[i].el] + cov[P[j].el])) bonds.push([i, j]);
  }

  var reduce = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  var rotY = -0.5, rotX = 0.35, dragging = false, lastX = 0, lastY = 0, idle = 0;
  var dpr = Math.min(window.devicePixelRatio || 1, 2);
  var W = 560, H = 520;
  canvas.width = W * dpr; canvas.height = H * dpr;

  function isDark() {
    var t = document.documentElement.getAttribute("data-theme");
    if (t) return t === "dark";
    return window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches;
  }

  function project() {
    var cy_ = Math.cos(rotY), sy = Math.sin(rotY), cx_ = Math.cos(rotX), sx = Math.sin(rotX);
    return P.map(function (p) {
      var x = p.x * cy_ + p.z * sy;
      var z = -p.x * sy + p.z * cy_;
      var y = p.y * cx_ - z * sx;
      z = p.y * sx + z * cx_;
      var persp = CAM / (CAM - z);
      return { x: x * persp, y: y * persp, z: z, s: persp, el: p.el };
    });
  }

  function draw() {
    var dark = isDark();
    var Q = project();
    var scale = 1.0 * (Math.min(W, H) / 2) / (maxR * CAM / (CAM - maxR));
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, W, H);
    ctx.save();
    ctx.translate(W / 2, H / 2);

    // hydrogen bonds: dashed
    ctx.setLineDash([5, 6]);
    ctx.lineWidth = 2;
    ctx.strokeStyle = dark ? "rgba(147,178,255,0.75)" : "rgba(42,82,189,0.7)";
    hbonds.forEach(function (b) {
      var a = Q[b[0]], c = Q[b[1]];
      ctx.beginPath();
      ctx.moveTo(a.x * scale, -a.y * scale);
      ctx.lineTo(c.x * scale, -c.y * scale);
      ctx.stroke();
    });
    ctx.setLineDash([]);

    // painter's algorithm: bonds and atoms sorted back to front
    var items = [];
    bonds.forEach(function (b) { items.push({ t: "b", z: (Q[b[0]].z + Q[b[1]].z) / 2 - 0.3, a: Q[b[0]], c: Q[b[1]] }); });
    Q.forEach(function (q) { items.push({ t: "a", z: q.z, q: q }); });
    items.sort(function (m, n) { return m.z - n.z; });

    items.forEach(function (it) {
      if (it.t === "b") {
        var w = 0.17 * scale * (it.a.s + it.c.s) / 2;
        ctx.lineCap = "round";
        ctx.lineWidth = w + 2.5;
        ctx.strokeStyle = dark ? "#0b1520" : "#3a4652";
        ctx.beginPath();
        ctx.moveTo(it.a.x * scale, -it.a.y * scale);
        ctx.lineTo(it.c.x * scale, -it.c.y * scale);
        ctx.stroke();
        ctx.lineWidth = w;
        ctx.strokeStyle = dark ? "#8796a6" : "#aab5bf";
        ctx.stroke();
      } else {
        var q = it.q, st = style[q.el];
        var r = st.r * scale * q.s;
        var x = q.x * scale, y = -q.y * scale;
        var g = ctx.createRadialGradient(x - r * 0.35, y - r * 0.4, r * 0.1, x, y, r);
        g.addColorStop(0, "#ffffff");
        g.addColorStop(0.28, st.c);
        g.addColorStop(1, shade(st.c, 0.55));
        ctx.beginPath();
        ctx.arc(x, y, r, 0, Math.PI * 2);
        ctx.fillStyle = g;
        ctx.fill();
        ctx.lineWidth = 1.5;
        ctx.strokeStyle = dark ? "#0b1520" : "#33404d";
        ctx.stroke();
      }
    });
    ctx.restore();
  }

  function shade(hex, f) {
    var n = parseInt(hex.slice(1), 16);
    var r = Math.round(((n >> 16) & 255) * f), g = Math.round(((n >> 8) & 255) * f), b = Math.round((n & 255) * f);
    return "rgb(" + r + "," + g + "," + b + ")";
  }

  var last = 0;
  function tick(t) {
    var dt = last ? Math.min(t - last, 50) : 16;
    last = t;
    if (!dragging && !reduce) {
      idle += dt;
      if (idle > 1200) rotY += dt * 0.00025;
    }
    draw();
    if (visible) requestAnimationFrame(tick);
    else last = 0;
  }

  // pause when off screen
  var visible = true;
  if ("IntersectionObserver" in window) {
    new IntersectionObserver(function (entries) {
      var was = visible;
      visible = entries[0].isIntersecting;
      if (visible && !was) requestAnimationFrame(tick);
    }).observe(canvas);
  }

  canvas.addEventListener("pointerdown", function (e) {
    dragging = true; lastX = e.clientX; lastY = e.clientY;
    canvas.setPointerCapture(e.pointerId);
  });
  canvas.addEventListener("pointermove", function (e) {
    if (!dragging) return;
    rotY += (e.clientX - lastX) * 0.01;
    rotX += (e.clientY - lastY) * 0.01;
    rotX = Math.max(-1.4, Math.min(1.4, rotX));
    lastX = e.clientX; lastY = e.clientY;
    if (reduce) draw();
  });
  function end() { dragging = false; idle = 0; }
  canvas.addEventListener("pointerup", end);
  canvas.addEventListener("pointercancel", end);

  if (window.matchMedia) {
    var mq = window.matchMedia("(prefers-color-scheme: dark)");
    if (mq.addEventListener) mq.addEventListener("change", draw);
  }

  draw();
  requestAnimationFrame(tick);
})();
