# create output directory
mkdir -p txt_output_files

# Uncomment for generating mesh with cylinder radius varying from [0.01,0.1]
for radius in $(seq 0.01 0.002 0.10)
do
     echo "-------------------------------------------------------------------------------------------------"
     echo "Cylinder radius =: $radius"
     python 2D-mesh-generator.py "$radius" &> "txt_output_files/mesh_gen_for_r_$radius.txt"
done