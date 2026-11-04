# BlurRelax
So many different ways to do a mush/relax/smooth

Slide Mode:
One of the most annoying things with trying to even out geometry is that the relax will either pin the un-even verts on the border, or it won't pin them at all and you lose your shape.
I added "Slide" mode for pinning that fixes this issue.
![Pinning Demo](references/BlurRelax_PinVsSlide.gif)

Volume Preservation:
Keep some percentage of the volume when smoothing. Uses the Taubin Smoothing algorithm

Reprojection:
After every iteration, snap the vertices back to the closest point on the original mesh (or a *subdivided* version of the original mesh).
This is more expensive than volume preservation, but more accurate. However, tight creases can cause chatter

Delta Mush:
It's delta mush

Reference Mush:
Only smooth out the difference between the input and the originalGeometry

Interpolation of iterations:
Iterations are linearly interpolated. So an iteration of 1.5 is a blend of 1 and 2 iterations.



