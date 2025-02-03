import streamlit as st
import plotly.graph_objects as go
import numpy as np
from scales_func import *

# Title
st.title("SCALES Exposure Time Calculator")

# Subtitle - Instrument Config
st.subheader("Instrument Config")

# Mode Selection
mode = st.radio("Mode", ["Imager-TODO", "IFS Low Res", "IFS Mid Res - TODO" ],index=1)

# Filter Selection
target_filter = st.selectbox("Filter", ["K", "L", "M", "CH4", "ICE","SED"],index=0)

# Subtitle - Source/Scene
st.subheader("Source/Scene")

# Source Type Selection
source_type = st.radio("Source Type", ["Point Source", "Extended Object -TODO"],index=0)

# Input Fields
ab_mag = st.number_input("AB Mag", min_value=0.0, step=0.1, format="%.1f",value=14.0)
itime = st.number_input("Integration time", min_value=0.0, step=1.0, format="%.1f",value=30.0)
nframes = st.number_input("no of frames", min_value=1, step=1,value=1)

filters = {
"K" : (1.95,2.45,2.185,200),
"L" : (2.90,4.15,3.46,80),
"M" : (4.5,5.2,5.83,140),
"CH4" : (3.1,3.5,3.3,250),
"ICE" : (2.0,4.0,2.82,50),
"SED" : (2.0,5.0,3.16,35),   
}

filter=target_filter
st.write("### Filter:	"+filter+ " ###")
st.write("### Minimum wavelength:	"+str(filters[filter][0])+ " ###")
st.write("### Maximum wavelength:	"+str(filters[filter][1])+ " ###")
st.write("### R:   "+ str(filters[filter][3])+ " ###")

snr,wave=scales_etc_snr(filter,itime,nframes,ab_mag)
if "snr" not in st.session_state:
	st.session_state.snr = snr
x = wave
y = snr[:,54,54]

fig = go.Figure()
fig.add_trace(go.Scatter(x=x, y=y, mode='lines', name='SNR Curve'))

fig.update_layout(title="SNR of peak spaxel ",
				  xaxis_title="Wavelength",
				  yaxis_title="SNR",
				  template="plotly_white")

st.plotly_chart(fig)

# Streamlit app
st.subheader("3D SNR Cube")

# Slider to select the slice index along the third axis
slice_index = st.slider("Select Channel", 0, snr.shape[0]-1, 0)
st.write('Wavelength :',str(round(wave[slice_index],2)))
# Get the selected slice
current_slice = snr[slice_index,:,:]

# Create the Plotly heatmap
fig = go.Figure(data=go.Heatmap(
    z=current_slice,
    colorscale="Viridis",
    hovertemplate="X: %{x}<br>Y: %{y}<br>Value: %{z}<extra></extra>"  # Display cursor value
))

# Update layout
fig.update_layout(
    title=f"Wavlength Slice {round(wave[slice_index],2)} of SNR Data Cube",
    xaxis_title="Spaxel",
    yaxis_title="Spaxel",
    width=500,
    height=500
)

# Display the Plotly figure in Streamlit
st.plotly_chart(fig)


