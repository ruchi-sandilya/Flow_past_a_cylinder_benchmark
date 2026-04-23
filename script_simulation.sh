
# Case 1
# Uncomment for simulating with cylinder radius between [0.01,0.05)
# Used for laminar flow Re = 20
Re=20
for radius in $(seq 0.01 0.002 0.048)
do
    echo "-------------------------------------------------------------------------------------------------"
    echo "Cylinder radius =: $radius"
    echo "Reynolds number =: $Re"
    python flow_simulation_generate_images.py "$radius" "$Re" \
    &> "txt_output_files/vel_out_r_${radius}_Re_${Re}.txt"
done

# Case 2
# Uncomment for simulating with cylinder radius between (0.06,0.1]
# Used for periodic flow Re = 120
Re=120
for radius in $(seq 0.062 0.002 0.1)
do
     echo "-------------------------------------------------------------------------------------------------"
     echo "Cylinder radius =: $radius"
     echo "Reynolds number =: $Re"
     python flow_simulation_generate_images.py "$radius" "$Re" \
    &> "txt_output_files/vel_out_r_${radius}_Re_${Re}.txt"
done

# Case 3
# Uncomment for simulating Reynolds Number between [20,40] (laminar flow)
radius=0.05
for ((Re=20; Re<=40; Re=Re+1)); do
     echo "-------------------------------------------------------------------------------------------------"
     echo "Cylinder radius =: $radius"
     echo "Reynolds number =: $Re"

     python flow_simulation_generate_images.py "$radius" "$Re" \
        &> "txt_output_files/vel_out_r_${radius}_Re_${Re}.txt"
done

# Case 4
# Uncomment for simulating Reynolds Number between [20,120]
radius=0.05
for ((Re=100; Re<=120; Re=Re+1)); do
     echo "-------------------------------------------------------------------------------------------------"
     echo "Cylinder radius =: $radius"
     echo "Reynolds number =: $Re"

     python flow_simulation_generate_images.py "$radius" "$Re" \
        &> "txt_output_files/vel_out_r_${radius}_Re_${Re}.txt"
done