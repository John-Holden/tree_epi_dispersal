#!/bin/bash

# Generate a simulation based on frames
rate1=1
rate2=1
#ffmpeg -r $rate1 -start_number 0 -i temp_frames/%04d.png -c:v libx264 -r $rate2 -pix_fmt yuv420p sim.mp4
#ffmpeg -framerate 1 -pattern_type glob -i temp_frames/'*.png' -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4
ffmpeg -framerate 1 -pattern_type glob -i temp_frames/'*.png' -c:v libx264 -pix_fmt yuv420p sim.mp4

#rm -rf temp_frames/0*  % remove temp frames
# rm -rf numpy_dat/0* % remove raw data
echo "animation complete:"
echo "REMOVE UNWANTED: raw data DATA & temp_frames!!!!"


