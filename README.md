# BlurRelax
Why does maya not have a smooth or relax deformer?  They've got brushes, sure, but no deformers.
This one is pretty fast. I've got it working with SIMD under the hood.

Also includes pinning and sliding along hard edges/borders.
One of the most annoying things with trying to even out geometry is that the relax will either pin the un-even verts on the border, or it won't pin them at all and you lose your shape.
I added "Slide" mode for pinning that fixes this issue.

![Pinning Demo](https://tbttfox.github.io/Images/BlurRelax_PinVsSlide.gif)

In future versions, I plan to add user defined pin sets.
