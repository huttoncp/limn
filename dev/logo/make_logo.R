## man/figures/logo.png from the sticker artwork in this folder.
##
## The artwork is a pointy-top hex on a white square, 2048 x 2048 px. Two
## steps make the logo:
##
## 1. The letter R is made brighter and easier to read. It's painted as a thin
##    blue outline whose strokes are nearly the colour of the glow behind them.
##    The outline is found by colour in a small window around the R. The strokes
##    are the largest region it encloses. Both are brightened in HSV, the
##    strokes filled so the letter reads as solid. Nothing else changes.
##
## 2. Everything outside the hex's ring is made transparent, the sticker is
##    cropped to the hex and scaled for the README and the pkgdown home page.
##    This follows illume's dev/logo/make_logo.R, with the background white
##    instead of black. Along each row, what lies before the first pixel darker
##    than white, or after the last, is outside. The ring's outer edge is
##    anti-aliased against the white, so it's "un-matted": a pixel that is
##    fraction a of the ring's colour and 1 - a of white gets alpha a and the
##    ring's colour.
##
## Run from the package root:  Rscript dev/logo/make_logo.R
## Needs png and magick; neither is a dependency of the package.

src   <- "dev/logo/limn-hex-source.jpg"
out   <- "man/figures/logo.png"
width <- 480L     # px; displayed at 139 px tall in the README, so ~4x for sharp screens
thr   <- 240      # 0-255: a lowest channel below this belongs to the sticker

A <- as.integer(magick::image_data(magick::image_read(src), "rgb")) / 255   # rows x cols x 3
nr <- dim(A)[1]; nc <- dim(A)[2]

## ---- 1. the R ---------------------------------------------------------------
wr <- 296:450; wc <- 955:1102              # a window around the R
W <- A[wr, wc, , drop = FALSE]
hr <- length(wr); hc <- length(wc)
## the outline: bright, saturated blue, not the network's cyan nodes
outline <- W[, , 3] > 0.80 & W[, , 2] < 0.72 & (W[, , 3] - W[, , 1]) > 0.55
grow1 <- function(M) {                     # close one-pixel gaps
  G <- M
  for (dr in -1:1) for (dc in -1:1) {
    S <- matrix(FALSE, nrow(M), ncol(M))
    rs <- max(1, 1 + dr):min(nrow(M), nrow(M) + dr); cs <- max(1, 1 + dc):min(ncol(M), ncol(M) + dc)
    S[rs, cs] <- M[rs - dr, cs - dc]
    G <- G | S
  }
  G
}
wall <- grow1(outline)
## label the regions the outline encloses (4-connected flood fill)
lab <- matrix(0L, hr, hc); k <- 0L
for (r0 in seq_len(hr)) for (c0 in seq_len(hc)) {
  if (wall[r0, c0] || lab[r0, c0]) next
  k <- k + 1L; stack <- c(r0, c0); lab[r0, c0] <- k
  while (length(stack)) {
    n <- length(stack); r <- stack[n - 1]; cc <- stack[n]; stack <- stack[-c(n - 1, n)]
    for (d in list(c(-1, 0), c(1, 0), c(0, -1), c(0, 1))) {
      r2 <- r + d[1]; c2 <- cc + d[2]
      if (r2 >= 1 && r2 <= hr && c2 >= 1 && c2 <= hc && !wall[r2, c2] && !lab[r2, c2]) {
        lab[r2, c2] <- k; stack <- c(stack, r2, c2)
      }
    }
  }
}
edge_labels <- unique(c(lab[1, ], lab[hr, ], lab[, 1], lab[, hc]))
sizes <- tabulate(lab[lab > 0], nbins = k)
enclosed <- setdiff(which(sizes > 0), edge_labels)
strokes <- lab == enclosed[which.max(sizes[enclosed])]    # between the outer and inner outlines
letter <- strokes | wall
cat(sprintf("R: outline %d px, strokes %d px, enclosed regions %d\n",
            sum(outline), sum(strokes), length(enclosed)))
box3 <- function(M) {                      # a 3 x 3 mean, for soft edges
  P <- matrix(0, nrow(M) + 2, ncol(M) + 2); P[2:(nrow(M) + 1), 2:(ncol(M) + 1)] <- M
  S <- 0
  for (dr in 0:2) for (dc in 0:2) S <- S + P[(1 + dr):(nrow(M) + dr), (1 + dc):(ncol(M) + dc)]
  S / 9
}
m_f <- as.vector(pmin(1, box3(box3(letter * 1)) * 1.4))
m_o <- as.vector(pmin(1, box3(outline * 1) * 1.6))
hsv <- grDevices::rgb2hsv(t(matrix(W, ncol = 3)), maxColorValue = 1)
fill_rgb <- t(grDevices::col2rgb(grDevices::hsv(hsv["h", ], hsv["s", ] * 0.95,
                                                 pmin(1, hsv["v", ] * 1.35)))) / 255
line_rgb <- t(grDevices::col2rgb(grDevices::hsv(hsv["h", ], hsv["s", ] * 0.75,
                                                 pmin(1, hsv["v", ] * 1.15)))) / 255
new <- matrix(W, ncol = 3) * (1 - m_f) + fill_rgb * m_f
new <- new * (1 - m_o) + line_rgb * m_o
A[wr, wc, ] <- array(new, c(hr, hc, 3))

## ---- 2. the cut-out -----------------------------------------------------------
mn <- 255 * pmin(A[, , 1], A[, , 2], A[, , 3])
lum <- 255 * (0.2126 * A[, , 1] + 0.7152 * A[, , 2] + 0.0722 * A[, , 3])

## OUTSIDE. The hex is convex, so along any row the sticker is one run of
## pixels: what lies before its first dark pixel or after its last is outside.
dark <- mn < thr
outside <- matrix(TRUE, nr, nc)
for (r in seq_len(nr)) {
  on <- which(dark[r, ])
  if (length(on)) outside[r, min(on):max(on)] <- FALSE
}

## chessboard distance from the edge, both ways, out to 3 px
shift <- function(M, dr, dc, fill) {
  out <- matrix(fill, nr, nc)
  rs <- max(1, 1 + dr):min(nr, nr + dr); cs <- max(1, 1 + dc):min(nc, nc + dc)
  out[rs, cs] <- M[rs - dr, cs - dc]
  out
}
grow <- function(M) {
  G <- M
  for (dr in -1:1) for (dc in -1:1) G <- G | shift(M, dr, dc, FALSE)
  G
}
d_in <- matrix(Inf, nr, nc)    # inside pixels: how far from the outside
d_out <- matrix(Inf, nr, nc)   # outside pixels: how far from the inside
Go <- outside; Gi <- !outside
for (k in 1:3) {
  Go <- grow(Go); d_in[Go & !outside & !is.finite(d_in)] <- k
  Gi <- grow(Gi); d_out[Gi & outside & !is.finite(d_out)] <- k
}

## the ring's own colour near each edge pixel: the mean of ring pixels 3 px
## in, within a small window. The ring's outer band is shaded -- about 140 at
## the top, 125 at the sides, 105 along the bottom -- so one colour won't do.
ring <- which(d_in == 3 & mn < 200)
band <- which((d_in <= 2) | (d_out <= 2))
rr <- ((band - 1) %% nr) + 1; cc <- ((band - 1) %/% nr) + 1
acc <- matrix(0, length(band), 3); cnt <- numeric(length(band))
is_ring <- logical(nr * nc); is_ring[ring] <- TRUE
for (dr in -4:4) for (dc in -4:4) {
  r2 <- rr + dr; c2 <- cc + dc
  ok <- r2 >= 1 & r2 <= nr & c2 >= 1 & c2 <= nc
  idx <- (c2[ok] - 1) * nr + r2[ok]
  hit <- is_ring[idx]
  w <- which(ok)[hit]
  for (ch in 1:3) acc[w, ch] <- acc[w, ch] + A[, , ch][idx[hit]]
  cnt[w] <- cnt[w] + 1
}
has <- cnt > 0
ringcol <- acc[has, , drop = FALSE] / cnt[has]
ringlum <- 255 * (0.2126 * ringcol[, 1] + 0.7152 * ringcol[, 2] + 0.0722 * ringcol[, 3])

## assemble RGBA: opaque inside, clear outside, un-matted along the edge
alpha <- matrix(0, nr, nc); alpha[!outside] <- 1
R <- A[, , 1]; G <- A[, , 2]; B <- A[, , 3]
bi <- band[has]
alpha[bi] <- pmin(pmax((255 - lum[bi]) / (255 - ringlum), 0), 1)
R[bi] <- ringcol[, 1]; G[bi] <- ringcol[, 2]; B[bi] <- ringcol[, 3]
cat(sprintf("edge pixels un-matted: %d (%d without a ring colour nearby, left as they were)\n",
            length(bi), sum(!has)))

## crop to the sticker, with a 2 px margin
keep <- which(alpha > 0.01, arr.ind = TRUE)
r0 <- max(1, min(keep[, 1]) - 2); r1 <- min(nr, max(keep[, 1]) + 2)
c0 <- max(1, min(keep[, 2]) - 2); c1 <- min(nc, max(keep[, 2]) + 2)
rgba <- array(0, c(r1 - r0 + 1, c1 - c0 + 1, 4))
rgba[, , 1] <- R[r0:r1, c0:c1]; rgba[, , 2] <- G[r0:r1, c0:c1]
rgba[, , 3] <- B[r0:r1, c0:c1]; rgba[, , 4] <- alpha[r0:r1, c0:c1]
cat(sprintf("cropped to %d x %d px (rows %d-%d, cols %d-%d)\n",
            c1 - c0 + 1, r1 - r0 + 1, r0, r1, c0, c1))

full <- tempfile(fileext = ".png")
png::writePNG(rgba, full)
img <- magick::image_read(full)
img <- magick::image_resize(img, sprintf("%dx", width), filter = "Lanczos")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
magick::image_write(img, out, format = "png", quality = 95)
info <- magick::image_info(magick::image_read(out))
cat(sprintf("wrote %s: %d x %d px, %.0f KB\n", out, info$width, info$height,
            file.size(out) / 1024))
