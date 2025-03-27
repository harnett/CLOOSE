import socket
import numpy as np
import h5py

# Path to your MATLAB v7.3 .mat file
mat_file = r"C:\Users\MysticSocialist\Documents\GitHub\CLOOSE_MS\Generate_images\WelcomeMsg.mat"

# Load the .mat file using h5py (MATLAB v7.3 files are HDF5-based)
with h5py.File(mat_file, 'r') as f:
    # Adjust the key if necessary. Here we assume the variable is named 'stack'.
    stack = np.array(f['stack'])

print("Original stack shape:", stack.shape)
# Reorder axes from (500, 796, 512) to (512, 796, 500)
stack = np.transpose(stack, (2, 1, 0))
print("Transposed stack shape:", stack.shape)

# Define TCP/IP communication parameters
host = "127.0.0.1"    # Loopback address (localhost)
port = 30001          # Port number (must match MATLAB receiver)

# Create a TCP server socket
with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as server_socket:
    server_socket.bind((host, port))
    server_socket.listen(1)
    print(f"Python TCP server listening on {host}:{port}...")
    
    # Wait for a connection from the MATLAB client
    conn, addr = server_socket.accept()
    with conn:
        print(f"Connected by {addr}")
        
        # Determine the number of frames (now expected shape: [512, 796, 500])
        num_frames = stack.shape[2]
        for iframe in range(num_frames):
            # Extract the 2D image frame (512 x 796)
            frame_2d = stack[:, :, iframe]
            
            # Flatten the frame in column-major (Fortran) order to match MATLAB's reshape behavior
            # Convert to little-endian 16-bit unsigned integers ('<u2')
            frame_flat = frame_2d.flatten(order='F').astype('<u2')
            
            # Send the raw bytes of the frame over the TCP connection
            conn.sendall(frame_flat.tobytes())
        
        print("All frames have been sent.")
