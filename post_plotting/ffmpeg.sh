#!/bin/bash

ffmpeg -framerate 24 -i output2/trace_%03d.png -vcodec libx264 -pix_fmt yuv420p -crf 10 ssva.mp4

# -crf: quality: the lower the better
