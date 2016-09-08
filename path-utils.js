// NEXT: post on red, integrate with p5.Font

/**
 * Gets an array of points along a text path
 *
 *  TODO: add params
 *
 * @param {object} options - an optional object that can contain:
 * <br>sampleFactor - the ratio of path-length to number of samples
 * (default=.25); higher values yield more points and are therefore
 * more precise
 * <br>simplifyThreshold - if set to a non-zero value, collinear points will be
 * be removed from the polygon; the value represents the threshold angle to use
 * when determining whether two edges are collinear
 * @return {array}  an array of points, each with x, y,
 * and alpha (the path angle)
 */
function textToPoints(font, txt, x, y, fontSize, options) {

  var xoff = 0, result = [], glyphs = font._getGlyphs(txt);

  for (var i = 0; i < glyphs.length; i++) {
    var gpath = glyphs[i].getPath(x, y, fontSize),
      paths = splitPaths(gpath.commands);

    for (var j = 0; j < paths.length; j++) {
      var pts = pathToPoints(paths[j], options);
      for (var k = 0; k < pts.length; k++) {
        pts[k].x += xoff;
        result.push(pts[k]);
      }
    }

    xoff += glyphs[i].advanceWidth * font._scale(fontSize);
  }

  return result;
}

/**
 * Gets an array of hierarchical Poly objects for the glyph
 * @param  {glyph} glyph the opentype glyph
 *
 * TODO: add params
 *
 * @param {object} options - an optional object that can contain:
 * <br>sampleFactor - the ratio of path-length to number of samples
 * (default=.25); higher values yield more points and are therefore
 * more precise
 * <br>simplifyThreshold - if set to a non-zero value, collinear points will be
 * be removed from the polygon; the value represents the threshold angle to use
 * when determining whether two edges are collinear
 * @return {array}  an array of hierarchical Poly objects
 */
function glyphToPolys(glyph, x, y, fontSize, options) {

  x = x || 0;
  y = y || 0;
  fontSize = fontSize || 100;

  var opts = parseOpts(options, {
    sampleFactor: .15,
    simplifyThreshold: 0,
  });

  var gpath = glyph.getPath(x, y, fontSize),
    paths = splitPaths(gpath.commands),
    polys = [];

  for (var j = 0; j < paths.length; j++) {

    var pts = pathToPoints(paths[j], opts);
    polys.push(new Poly(pts));
  }

  for (var j = 0; j < polys.length; j++) {

    // check each of the other polys to see if we are a hole
    for (var i = 0; i < polys.length; i++) {

      //console.log('Checking '+j+' and '+i);
      if (i != j && polyInPoly(polys[i], polys[j])) {

        polys[i].addHole(polys[j].points);
        polys.splice(j--, 1);
        break;
      }
    }
  }

  return polys;
}

function pathToPoints(cmds, options) {

  var opts = parseOpts(options, {
    sampleFactor: .1,
    simplifyThreshold: 0,
  });

  var len = pointAtLength(cmds,0,1), // total-length
    t = len / (len * opts.sampleFactor),
    pts = [];

  for (var i = 0; i < len; i += t)
    pts.push(pointAtLength(cmds, i));

  if (opts.simplifyThreshold) {
    var count = simplify(pts, opts.simplifyThreshold);
    //console.log('Simplify: removed ' + count + ' pts');
  }

  return pts;
}

function simplify(pts, angle) {

  if (typeof angle === 'undefined') angle = 0;

  var num = 0;
  for (var i = pts.length - 1; pts.length > 3 && i >= 0; --i) {

    if (collinear(at(pts, i - 1), at(pts, i), at(pts, i + 1), angle)) {

      // Remove the middle point
      pts.splice(i % pts.length, 1);
      num++;
    }
  }
  return num;
}

function pathBounds(pts) {

  var minX = minY = Number.MAX_VALUE;
  var maxX = maxY = -Number.MAX_VALUE;

  for (var i = 0; i < pts.length; i++) {

    if (pts[i].x < minX) minX = pts[i].x;
    if (pts[i].x > maxX) maxX = pts[i].x;

    if (pts[i].y < minY) minY = pts[i].y;
    if (pts[i].y > maxY) maxY = pts[i].y;
  }

  var w = maxX - minX,
    h = maxY - minY;

  return {
    x: minX, y: minY, x2: maxX, y2: maxY,
    w: w, h: h, cx: minX + w / 2, cy: minY + h / 2
  };
}

function splitPaths(cmds) {

  var paths = [], current;
  for (var i = 0; i < cmds.length; i++) {
    if (cmds[i].type === 'M') {
      current && paths.push(current);
      current = [];
    }
    current.push(cmdToArr(cmds[i]));
  }
  paths.push(current);

  return paths;
}

function cmdToArr(cmd) {

  var arr = [ cmd.type ];
  if (cmd.type === 'M' || cmd.type === 'L') { // moveto or lineto
    arr.push(cmd.x, cmd.y);
  } else if (cmd.type === 'C') {
    arr.push(cmd.x1, cmd.y1, cmd.x2, cmd.y2, cmd.x, cmd.y);
  } else if (cmd.type === 'Q') {
    arr.push(cmd.x1, cmd.y1, cmd.x, cmd.y);
  }
  // else if (cmd.type === 'Z') { /* no-op */ }
  return arr;
}

/*
 * Checks if point is within a (simple) polygon, as specified by its vertices
 * (note: does not handle holes)
 */
function pointInPoly(x, y, pts, bbox) {

  // ray-casting based on
  // http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html

  var inside = false;
  if (!bbox || pointInBox(x, y, bbox)) {

    for (var i = 0, j = pts.length - 1; i < pts.length; j = i++) {

      var xi = pts[i].x,
        yi = pts[i].y,
        xj = pts[j].x,
        yj = pts[j].y;

      var intersect = ((yi > y) != (yj > y)) && (x < (xj - xi) *
        (y - yi) / (yj - yi) + xi);

      if (intersect) inside = !inside;
    }
  }
  return inside;
}

/*
 * Naive check if polygon2 is inside polygon1 (should also check
 * 	intersection and holes)
 * @param polygon1 polygon that contains other
 * @param polygon2 polygon that may be inner to other
 * @return true if polygon2 is inside polygon1, else false
 */
function polyInPoly(p1, p2) {

  if (!Array.isArray(p1) && typeof p1.points !== 'undefined')
    p1 = p1.points;

  if (!Array.isArray(p2) && typeof p2.points !== 'undefined')
    p2 = p2.points;

  var bb1 = pathBounds(p1),
    bb2 = pathBounds(p2),
    pIn = pointInBox;

  // first check 4 corners of bbox2 of are inside bbox1
  if (!(pIn(bb2.x, bb2.y, bb1) && pIn(bb2.x2, bb2.y, bb1) &&
      pIn(bb2.x, bb2.y2, bb1) && pIn(bb2.x2, bb2.y2, bb1))) {
    return false;
  }

  // if true, then check each point is inside
  for (var i = 0; i < p2.length; i++) {
    if (!pointInPoly(p2[i].x, p2[i].y, p1))
      return false;
  }

  return true;
}

function pointInBox(x, y, bbox) {

  return x >= bbox.x && x <= bbox.x2 && y >= bbox.y && y <= bbox.y2;
}

function joinBounds(bbs) {

  var minX = minY = Number.MAX_VALUE,
    maxX = maxY = -Number.MAX_VALUE;

  for (var i = 0; i < bbs.length; i++) {

    if (bbs[i].x < minX) minX = bbs[i].x;
    if (bbs[i].x2 > maxX) maxX = bbs[i].x2;

    if (bbs[i].y < minY) minY = bbs[i].y;
    if (bbs[i].y2 > maxY) maxY = bbs[i].y2;
  }

  return {
    x: minX, y: minY, w: maxX - minX, h: maxY - minY, x2: maxX, y2: maxY
  }
}

function parseOpts(options, defaults) {

  if (typeof options !== 'object') {
    options = defaults;
  }
  else {
    for (var key in defaults) {
      if (typeof options[key] === 'undefined')
        options[key] = defaults[key];
    }
  }
  return options;
}

//////////////////////// class Poly //////////////////////////

function Poly(pts) {

  this.points = pts; // 1d array of contour point objects
  this.holes = undefined; // 2d array of hole point objects
}

Poly.prototype.addHole = function(pts) {

  this.holes = this.holes || [];
  this.holes.push(pts);
}

Poly.prototype.bounds = function() {

  return pathBounds(this.points);
}

Poly.prototype.simplify = function(angle) {

  return simplify(this.points, angle);
};

/**
 * Tests for containment, either of a point, or another Poly
 * @param a point or another Poly object
 * @return {boolean}
 */
Poly.prototype.contains = function() {

  if (arguments.length > 1 && typeof arguments[0].x != 'undefined')
    return this.containsPoint(arguments[0], arguments[1]);

  return polyInPoly(this, arguments[0]);
}

Poly.prototype.containsPoint = function(x, y) {

  if (!pointInBox(x, y, pathBounds(this.points)))
    return false;

  if (!pointInPoly(x, y, this.points))
    return false;

  if (this.holes) {
    for (var i = 0; i < this.holes.length; i++) {
      if (pointInPoly(x, y, this.holes[i]))
        return false;
    }
  }

  return true;
}

Poly.prototype.getPoints = function(options) {

  var pts = [],
    opts = parseOpts(options, { includeHoles: true });

  for (var j = 0; j < this.points.length; j++) {
    pts.push(this.points[j]);
  }

  if (this.holes && opts.includeHoles) { // TODO: test/fix

    //console.log('poly returning holes');
    for (var i = 0; i < this.holes.length; i++) {
      for (var k = 0; k < this.holes[i].length; k++)
        pts.push(this.holes[i][k]);
    }
  }
  return pts;
}

Poly.prototype.draw = function() { // requires p5js

  beginShape();
  for (var j = 0; j < this.points.length; j++) {
    vertex(this.points[j].x, this.points[j].y);
  }
  if (this.holes) {
    for (var i = 0; i < this.holes.length; i++) {
      beginContour();
      var pts = this.holes[i];
      for (var k = 0; k < pts.length; k++)
        vertex(pts[k].x, pts[k].y);
      endContour();
    }
  }
  endShape(CLOSE);
}

Poly.prototype.drawPoints = function(pointSz) { // requires p5js.ellipse

  var s = pointSz || 2, pts = this.getPoints();
  for (var j = 0; j < pts.length; j++) {
    ellipse(pts[j].x, pts.y, s, s);
  }
}

//////////////////////// Helpers ////////////////////////////

function at(v, i) {
  var s = v.length;
  return v[i < 0 ? i % s + s : i % s];
}

function collinear(a, b, c, thresholdAngle) {

  if (!thresholdAngle) return areaTriangle(a, b, c) == 0;

  if (typeof collinear.tmpPoint1 === 'undefined') {
    collinear.tmpPoint1 = [];
    collinear.tmpPoint2 = [];
  }

  var ab = collinear.tmpPoint1, bc = collinear.tmpPoint2;
  ab.x = b.x - a.x;
  ab.y = b.y - a.y;
  bc.x = c.x - b.x;
  bc.y = c.y - b.y;

  var dot = ab.x * bc.x + ab.y * bc.y,
    magA = Math.sqrt(ab.x * ab.x + ab.y * ab.y),
    magB = Math.sqrt(bc.x * bc.x + bc.y * bc.y),
    angle = Math.acos(dot / (magA * magB));

  return angle < thresholdAngle;
};

function areaTriangle(a, b, c) {
  return (((b[0] - a[0]) * (c[1] - a[1])) - ((c[0] - a[0]) * (b[1] - a[1])));
}

/*
Portions of the code below are copyright © 2008 Dmitry Baranovskiy (via MIT license)
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

function findDotsAtSegment(p1x, p1y, c1x, c1y, c2x, c2y, p2x, p2y, t) {

  var t1 = 1 - t, t13 = pow(t1, 3), t12 = pow(t1, 2), t2 = t * t, t3 = t2 * t,
    x = t13 * p1x + t12 * 3 * t * c1x + t1 * 3 * t * t * c2x + t3 * p2x,
    y = t13 * p1y + t12 * 3 * t * c1y + t1 * 3 * t * t * c2y + t3 * p2y,
    mx = p1x + 2 * t * (c1x - p1x) + t2 * (c2x - 2 * c1x + p1x),
    my = p1y + 2 * t * (c1y - p1y) + t2 * (c2y - 2 * c1y + p1y),
    nx = c1x + 2 * t * (c2x - c1x) + t2 * (p2x - 2 * c2x + c1x),
    ny = c1y + 2 * t * (c2y - c1y) + t2 * (p2y - 2 * c2y + c1y),
    ax = t1 * p1x + t * c1x, ay = t1 * p1y + t * c1y,
    cx = t1 * c2x + t * p2x, cy = t1 * c2y + t * p2y,
    alpha = (90 - Math.atan2(mx - nx, my - ny) * 180 / PI);

  (mx > nx || my < ny) && (alpha += 180);

  return {
    x: x, y: y, m: { x: mx, y: my }, n: { x: nx, y: ny },
    start: { x: ax, y: ay }, end: { x: cx, y: cy }, alpha: alpha
  };
}

function getPointAtSegmentLength(p1x, p1y,
  c1x, c1y, c2x, c2y, p2x, p2y, length) {
  return (length == null) ? bezlen(p1x, p1y, c1x, c1y, c2x, c2y, p2x, p2y) :
    findDotsAtSegment(p1x, p1y, c1x, c1y, c2x, c2y, p2x, p2y,
      getTatLen(p1x, p1y, c1x, c1y, c2x, c2y, p2x, p2y, length));
}

function pointAtLength(path, length, istotal) {

  path = path2curve(path);
  var x, y, p, l, sp = '', subpaths = {}, point, len = 0;
  for (var i = 0, ii = path.length; i < ii; i++) {
    p = path[i];
    if (p[0] == 'M') {
      x = +p[1];
      y = +p[2];
    } else {
      l = getPointAtSegmentLength(x, y, p[1], p[2], p[3], p[4], p[5], p[6]);
      if (len + l > length) {
        if (!istotal) {
          point = getPointAtSegmentLength(x, y, p[1], p[2], p[3], p[4], p[5],
            p[6], length - len);
          return { x: point.x, y: point.y, alpha: point.alpha };
        }
      }
      len += l;
      x = +p[5];
      y = +p[6];
    }
    sp += p.shift() + p;
  }
  subpaths.end = sp;
  point = istotal ? len : subpath ? subpaths : findDotsAtSegment(x, y,
    p[0], p[1], p[2], p[3], p[4], p[5], 1);
  point.alpha && (point = { x: point.x, y: point.y, alpha: point.alpha });
  return point;
}

function pathToAbsolute(pathArray) {

  var res = [], x = 0, y = 0, mx = 0, my = 0, start = 0;
  if (pathArray[0][0] == 'M') {
    x = +pathArray[0][1];
    y = +pathArray[0][2];
    mx = x;
    my = y;
    start++;
    res[0] = ['M', x, y];
  }

  var crz = pathArray.length == 3 && pathArray[0][0] == 'M' &&
    pathArray[1][0].toUpperCase() == 'R' &&
    pathArray[2][0].toUpperCase() == 'Z';

  for (var r, pa, i = start, ii = pathArray.length; i < ii; i++) {
    res.push(r = []);
    pa = pathArray[i];
    if (pa[0] != String.prototype.toUpperCase.call(pa[0])) {
      r[0] = String.prototype.toUpperCase.call(pa[0]);
      switch (r[0]) {
        case 'A':
          r[1] = pa[1];
          r[2] = pa[2];
          r[3] = pa[3];
          r[4] = pa[4];
          r[5] = pa[5];
          r[6] = +(pa[6] + x);
          r[7] = +(pa[7] + y);
          break;
        case 'V':
          r[1] = +pa[1] + y;
          break;
        case 'H':
          r[1] = +pa[1] + x;
          break;
        case 'R':
          var dots = [x, y].concat(pa.slice(1));
          for (var j = 2, jj = dots.length; j < jj; j++) {
            dots[j] = +dots[j] + x;
            dots[++j] = +dots[j] + y;
          }
          res.pop();
          res = res.concat(catmullRom2bezier(dots, crz));
          break;
        case 'M':
          mx = +pa[1] + x;
          my = +pa[2] + y;
        default:
          for (j = 1, jj = pa.length; j < jj; j++) {
            r[j] = +pa[j] + ((j % 2) ? x : y);
          }
      }
    } else if (pa[0] == 'R') {
      dots = [x, y].concat(pa.slice(1));
      res.pop();
      res = res.concat(catmullRom2bezier(dots, crz));
      r = ['R'].concat(pa.slice(-2));
    } else {
      for (var k = 0, kk = pa.length; k < kk; k++) {
        r[k] = pa[k];
      }
    }
    switch (r[0]) {
      case 'Z':
        x = mx;
        y = my;
        break;
      case 'H':
        x = r[1];
        break;
      case 'V':
        y = r[1];
        break;
      case 'M':
        mx = r[r.length - 2];
        my = r[r.length - 1];
      default:
        x = r[r.length - 2];
        y = r[r.length - 1];
    }
  }

  return res;
}

function path2curve(path, path2) {

  var p = pathToAbsolute(path), p2 = path2 && pathToAbsolute(path2),
    attrs = { x: 0, y: 0, bx: 0, by: 0, X: 0, Y: 0, qx: null, qy: null },
    attrs2 = { x: 0, y: 0, bx: 0, by: 0, X: 0, Y: 0, qx: null, qy: null },

    processPath = function(path, d, pcom) {
      var nx, ny, tq = { T: 1, Q: 1 };
      if (!path) {
        return ['C', d.x, d.y, d.x, d.y, d.x, d.y];
      }!(path[0] in tq) && (d.qx = d.qy = null);
      switch (path[0]) {
        case 'M':
          d.X = path[1];
          d.Y = path[2];
          break;
        case 'A':
          path = ['C'].concat(a2c[apply](0, [d.x, d.y].concat(path.slice(1))));
          break;
        case 'S':
          if (pcom == 'C' || pcom == 'S') {
            nx = d.x * 2 - d.bx;
            ny = d.y * 2 - d.by;
          } else {
            nx = d.x;
            ny = d.y;
          }
          path = ['C', nx, ny].concat(path.slice(1));
          break;
        case 'T':
          if (pcom == 'Q' || pcom == 'T') {
            d.qx = d.x * 2 - d.qx;
            d.qy = d.y * 2 - d.qy;
          } else {
            d.qx = d.x;
            d.qy = d.y;
          }
          path = ['C'].concat(q2c(d.x, d.y, d.qx, d.qy, path[1], path[2]));
          break;
        case 'Q':
          d.qx = path[1];
          d.qy = path[2];
          path = ['C'].concat(q2c(d.x, d.y, path[1], path[2], path[3], path[4]));
          break;
        case 'L':
          path = ['C'].concat(l2c(d.x, d.y, path[1], path[2]));
          break;
        case 'H':
          path = ['C'].concat(l2c(d.x, d.y, path[1], d.y));
          break;
        case 'V':
          path = ['C'].concat(l2c(d.x, d.y, d.x, path[1]));
          break;
        case 'Z':
          path = ['C'].concat(l2c(d.x, d.y, d.X, d.Y));
          break;
      }
      return path;
    },

    fixArc = function(pp, i) {
      if (pp[i].length > 7) {
        pp[i].shift();
        var pi = pp[i];
        while (pi.length) {
          pcoms1[i] = 'A';
          p2 && (pcoms2[i] = 'A');
          pp.splice(i++, 0, ['C'].concat(pi.splice(0, 6)));
        }
        pp.splice(i, 1);
        ii = Math.max(p.length, p2 && p2.length || 0);
      }
    },

    fixM = function(path1, path2, a1, a2, i) {
      if (path1 && path2 && path1[i][0] == 'M' && path2[i][0] != 'M') {
        path2.splice(i, 0, ['M', a2.x, a2.y]);
        a1.bx = 0;
        a1.by = 0;
        a1.x = path1[i][1];
        a1.y = path1[i][2];
        ii = Math.max(p.length, p2 && p2.length || 0);
      }
    },

    pcoms1 = [], // path commands of original path p
    pcoms2 = [], // path commands of original path p2
    pfirst = '', // temporary holder for original path command
    pcom = ''; // holder for previous path command of original path

  for (var i = 0, ii = Math.max(p.length, p2 && p2.length || 0); i < ii; i++) {
    p[i] && (pfirst = p[i][0]); // save current path command

    if (pfirst != 'C') {
      pcoms1[i] = pfirst; // Save current path command
      i && (pcom = pcoms1[i - 1]); // Get previous path command pcom
    }
    p[i] = processPath(p[i], attrs, pcom);

    if (pcoms1[i] != 'A' && pfirst == 'C') pcoms1[i] = 'C';

    fixArc(p, i); // fixArc adds also the right amount of A:s to pcoms1

    if (p2) { // the same procedures is done to p2
      p2[i] && (pfirst = p2[i][0]);
      if (pfirst != 'C') {
        pcoms2[i] = pfirst;
        i && (pcom = pcoms2[i - 1]);
      }
      p2[i] = processPath(p2[i], attrs2, pcom);

      if (pcoms2[i] != 'A' && pfirst == 'C') pcoms2[i] = 'C';

      fixArc(p2, i);
    }
    fixM(p, p2, attrs, attrs2, i);
    fixM(p2, p, attrs2, attrs, i);
    var seg = p[i], seg2 = p2 && p2[i], seglen = seg.length,
      seg2len = p2 && seg2.length;
    attrs.x = seg[seglen - 2];
    attrs.y = seg[seglen - 1];
    attrs.bx = parseFloat(seg[seglen - 4]) || attrs.x;
    attrs.by = parseFloat(seg[seglen - 3]) || attrs.y;
    attrs2.bx = p2 && (parseFloat(seg2[seg2len - 4]) || attrs2.x);
    attrs2.by = p2 && (parseFloat(seg2[seg2len - 3]) || attrs2.y);
    attrs2.x = p2 && seg2[seg2len - 2];
    attrs2.y = p2 && seg2[seg2len - 1];
  }

  return p2 ? [p, p2] : p;
}

// http://schepers.cc/getting-to-the-point
function catmullRom2bezier(crp, z) {
  var d = [];
  for (var i = 0, iLen = crp.length; iLen - 2 * !z > i; i += 2) {
    var p = [{
      x: +crp[i - 2],
      y: +crp[i - 1]
    }, {
      x: +crp[i],
      y: +crp[i + 1]
    }, {
      x: +crp[i + 2],
      y: +crp[i + 3]
    }, {
      x: +crp[i + 4],
      y: +crp[i + 5]
    }];
    if (z) {
      if (!i) {
        p[0] = {
          x: +crp[iLen - 2],
          y: +crp[iLen - 1]
        };
      } else if (iLen - 4 == i) {
        p[3] = {
          x: +crp[0],
          y: +crp[1]
        };
      } else if (iLen - 2 == i) {
        p[2] = {
          x: +crp[0],
          y: +crp[1]
        };
        p[3] = {
          x: +crp[2],
          y: +crp[3]
        };
      }
    } else {
      if (iLen - 4 == i) {
        p[3] = p[2];
      } else if (!i) {
        p[0] = {
          x: +crp[i],
          y: +crp[i + 1]
        };
      }
    }
    d.push(['C', (-p[0].x + 6 * p[1].x + p[2].x) / 6, (-p[0].y + 6 * p[1].y +
      p[2].y) / 6, (p[1].x + 6 * p[2].x - p[3].x) / 6, (p[1].y + 6 * p[2].y -
      p[3].y) / 6, p[2].x, p[2].y ]);
  }

  return d;
}

function l2c(x1, y1, x2, y2) { return [x1, y1, x2, y2, x2, y2]; }

function q2c(x1, y1, ax, ay, x2, y2) {
  var _13 = 1 / 3, _23 = 2 / 3;
  return [
    _13 * x1 + _23 * ax, _13 * y1 + _23 * ay,
    _13 * x2 + _23 * ax, _13 * y2 + _23 * ay, x2, y2
  ];
}

function bezlen(x1, y1, x2, y2, x3, y3, x4, y4, z) {
  if (z == null) z = 1;
  z = z > 1 ? 1 : z < 0 ? 0 : z;
  var z2 = z / 2,
    n = 12, Tvalues = [-0.1252, 0.1252, -0.3678, 0.3678, -0.5873, 0.5873, -0.7699, 0.7699, -0.9041, 0.9041, -0.9816, 0.9816],
    sum = 0, Cvalues = [0.2491, 0.2491, 0.2335, 0.2335, 0.2032, 0.2032, 0.1601, 0.1601, 0.1069, 0.1069, 0.0472, 0.0472 ];
  for (var i = 0; i < n; i++) {
    var ct = z2 * Tvalues[i] + z2,
      xbase = base3(ct, x1, x2, x3, x4),
      ybase = base3(ct, y1, y2, y3, y4),
      comb = xbase * xbase + ybase * ybase;
    sum += Cvalues[i] * Math.sqrt(comb);
  }
  return z2 * sum;
}

function getTatLen(x1, y1, x2, y2, x3, y3, x4, y4, ll) {
  if (ll < 0 || bezlen(x1, y1, x2, y2, x3, y3, x4, y4) < ll) return;
  var t = 1, step = t / 2, t2 = t - step, l, e = .01;
  l = bezlen(x1, y1, x2, y2, x3, y3, x4, y4, t2);
  while (abs(l - ll) > e) {
    step /= 2;
    t2 += (l < ll ? 1 : -1) * step;
    l = bezlen(x1, y1, x2, y2, x3, y3, x4, y4, t2);
  }
  return t2;
}

function base3(t, p1, p2, p3, p4) {
  var t1 = -3 * p1 + 9 * p2 - 9 * p3 + 3 * p4,
    t2 = t * t1 + 6 * p1 - 12 * p2 + 6 * p3;
  return t * t2 - 3 * p1 + 3 * p2;
}
