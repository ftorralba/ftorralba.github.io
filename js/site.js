// Shared sidebar and disclaimer. To add, rename or reorder a section, edit PAGES below;
// every page picks up the change. Each page only needs <nav id="nav"></nav> and
// <footer id="disclaimer"></footer> plus this script.

var PAGES = [
  { href: "index.html",         label: "About" },
  { href: "resume.html",        label: "Resume" },
  { href: "current-work.html",  label: "Current work" },
  { href: "publications.html",  label: "Academic publications" },
  { href: "other-publications.html", label: "Other publications" },
  { href: "working-papers.html", label: "Working papers" },
  { href: "code-and-data.html", label: "Code and data" }
];

var NAME = "Francisco Torralba";
var DISCLAIMER = "The material on this website does not represent the views of the " +
  "Federal Reserve Bank of New York or the Federal Reserve System.";

(function () {
  var current = location.pathname.split("/").pop() || "index.html";

  var nav = document.getElementById("nav");
  if (nav) {
    var html = '<a class="site-name" href="index.html">' + NAME + '</a><ul>';
    PAGES.forEach(function (p) {
      var cur = p.href === current ? ' aria-current="page"' : "";
      html += '<li><a href="' + p.href + '"' + cur + ">" + p.label + "</a></li>";
    });
    nav.innerHTML = html + "</ul>";
  }

  var foot = document.getElementById("disclaimer");
  if (foot) foot.textContent = DISCLAIMER;
})();
