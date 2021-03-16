#!/bin/bash

# Generate a simulation based on frames
rate1=20
rate2=20
ls temp_frames

ffmpeg -r 5 -i temp_frames/%04d.png -c:v libx264 -vf fps=25 -pix_fmt yuv420p out.mp4

rm -rf temp_frames/*.png
echo "animation complete:"
echo "REMOVE UNWANTED: raw data DATA & temp_frames!!!!"


