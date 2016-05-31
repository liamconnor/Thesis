#!/bin/bash

montage -mode concatenate -tile x1 dm/cproj_p2_0000.jpg nu/cproj_p2_0000_nu.jpg final/slice.jpg 
composite -compose atop -dissolve 100 -geometry +80+50 -gravity northwest labels/z0.0.pdf final/slice.jpg final/slice.jpg
composite -compose atop -dissolve 100 -geometry +90+70 -gravity southwest labels/hline.pdf final/slice.jpg final/slice.jpg
composite -compose atop -dissolve 100 -geometry +95+90 -gravity southwest labels/length.pdf final/slice.jpg final/slice.jpg 

