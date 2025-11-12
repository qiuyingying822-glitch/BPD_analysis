
#!/bin/bash

# SparCC Installation Script
# This script installs SparCC and its dependencies

echo "=== Installing SparCC ==="

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo "Error: Python3 is required but not installed"
    exit 1
fi

# Create installation directory
INSTALL_DIR="$HOME/sparcc_installation"
mkdir -p "$INSTALL_DIR"
cd "$INSTALL_DIR"

echo "1. Downloading SparCC..."
# Clone SparCC repository (if available)
if command -v git &> /dev/null; then
    git clone https://github.com/attedgi/SparCC.git sparcc_repo
    SPARCC_DIR="$INSTALL_DIR/sparcc_repo"
else
    echo "Warning: git not available, please download SparCC manually"
    echo "Download from: https://bitbucket.org/yonatanf/sparcc"
    echo "And extract to: $INSTALL_DIR"
    exit 1
fi

echo "2. Setting up environment..."
# Add SparCC to PATH
echo "export PATH=\"$SPARCC_DIR:\$PATH\"" >> "$HOME/.bashrc"
export PATH="$SPARCC_DIR:$PATH"

echo "3. Testing installation..."
if command -v SparCC.py &> /dev/null; then
    echo "SparCC installation successful!"
    echo "SparCC location: $SPARCC_DIR/SparCC.py"
else
    echo "Error: SparCC installation failed"
    echo "Please check the installation manually"
    exit 1
fi

echo "4. Installing Python dependencies..."
# Install required Python packages
pip3 install numpy scipy

echo "=== SparCC Installation Completed ==="
echo ""
echo "To use SparCC, either:"
echo "1. Restart your terminal, or"
echo "2. Run: source ~/.bashrc"
echo ""
echo "Test the installation with: SparCC.py --help"
