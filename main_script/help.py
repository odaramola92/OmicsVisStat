#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OmicsVisStat - Help System
Standalone help script that can be run independently or imported by the GUI.

Usage:
    python help.py                    # Show help menu
    python help.py --tab statistics   # Show specific tab help
    from help import get_help_content # Import in GUI
"""

import sys
import io

# Set UTF-8 encoding for console output
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

def get_help_content():
    """
    Returns comprehensive help content for all tabs in OmicsVisStat.
    
    Returns:
        dict: Dictionary with tab names as keys and help content as values
    """
    return {
        "Statistics": """📊 STATISTICS TAB

Perform statistical analysis on metabolite data.

OVERVIEW:
The Statistics tab provides comprehensive tools for analyzing metabolite abundance, fold changes, and statistical significance.

KEY FEATURES:
• Differential expression analysis
• Volcano plots
• PCA (Principal Component Analysis)
• Heatmap generation
• Statistical summaries and reporting
• Multiple testing correction
• Data normalization options

HOW TO USE:
1. Load your metabolite data (Excel file)
2. Select analysis type:
   - Differential Expression
   - PCA
   - Clustering
3. Configure statistical parameters:
   - Significance threshold (p-value)
   - Fold change cutoff
   - Correction method
4. Run analysis
5. Review results and plots
6. Export results and visualizations

STATISTICAL METHODS:
• T-test (paired/unpaired)
• ANOVA (one-way, two-way)
• Multiple testing correction:
  - False Discovery Rate (FDR/Benjamini-Hochberg)
  - Bonferroni correction
  - Holm-Bonferroni
• Fold change calculations (linear or log2)

FILTERING OPTIONS:
• P-value threshold: 0.05, 0.01, 0.001, or custom
• Log2FC threshold: ±0.5, ±1.0, ±1.5, or custom
• Minimum detection threshold
• Sample size requirements

OUTPUTS:
• Statistical summary tables
• Significant metabolites list
• Volcano plot (log2FC vs -log10(p-value))
• PCA score plots
• Loading plots
• Heatmaps with hierarchical clustering

DATA REQUIREMENTS:
• Minimum 3 replicates per group (recommended)
• Numerical abundance/intensity values
• Missing value handling options
• Normalization recommended

TIPS:
⭐ Ensure sufficient sample size for statistical power
⭐ Apply appropriate multiple testing correction (FDR recommended)
⭐ Check data distribution before analysis (normality tests)
⭐ Use log-transformation for abundance data
⭐ Visualize data before and after normalization
⭐ Document analysis parameters for reproducibility
        """,
        
        "Visualization": """📈 VISUALIZATION TAB

Create publication-quality plots and charts.

OVERVIEW:
The Visualization tab generates various plots to visualize metabolite data and analysis results with professional styling and export options.

KEY FEATURES:
• Interactive plots with zoom/pan
• Multiple chart types
• Customizable styling and themes
• Export options (PNG, PDF, SVG, EPS)
• Color scheme selection
• Font and size customization
• Multi-plot layouts

AVAILABLE PLOTS:
• Volcano plots (differential expression)
• Bar charts (metabolite comparisons)
• Heatmaps (hierarchical clustering)
• PCA plots (2D and 3D)
• Box plots (distribution comparisons)
• Scatter plots (correlation analysis)
• Venn diagrams (overlap analysis)

HOW TO USE:
1. Load analyzed data
2. Select plot type from dropdown
3. Configure plot parameters:
   - X/Y axis variables
   - Color coding scheme
   - Point/bar sizes
   - Labels and titles
4. Preview the visualization
5. Adjust styling as needed
6. Export or save plot

CUSTOMIZATION OPTIONS:
• Color schemes: Viridis, Plasma, Set1, Set2, Custom
• Point markers: Circle, Square, Triangle, Diamond
• Line styles: Solid, Dashed, Dotted
• Font families: Arial, Times, Helvetica, Courier
• Font sizes: Title, axis labels, tick labels
• Figure dimensions: Width, height, DPI

EXPORT FORMATS:
• PNG: High-resolution raster (300+ DPI)
• PDF: Vector format for publications
• SVG: Scalable vector graphics (web/editing)
• EPS: Encapsulated PostScript (legacy journals)

INTERACTIVE FEATURES:
• Zoom: Box select to zoom in
• Pan: Drag to move view
• Reset: Double-click to reset view
• Hover: Show data point information
• Select: Click points to highlight

TIPS:
⭐ Use consistent color schemes across all plots
⭐ Add clear axis labels with units
⭐ Include informative titles
⭐ Export in vector format (PDF/SVG) for publications
⭐ Use 300+ DPI for high-quality PNG images
⭐ Save plot configurations for reproducibility
⭐ Test colorblind-friendly palettes
        """
    }


def display_help_menu():
    """Display interactive help menu in console."""
    help_content = get_help_content()
    
    print("\n" + "="*70)
    print("  OMICSVISSTAT - HELP SYSTEM")
    print("="*70 + "\n")
    
    print("Available Help Topics:\n")
    topics = list(help_content.keys())
    for i, topic in enumerate(topics, 1):
        print(f"  {i}. {topic}")
    
    print("\n  0. Exit")
    print("\n" + "-"*70)
    
    while True:
        try:
            choice = input("\nSelect a topic (0-{}): ".format(len(topics)))
            
            if choice == '0':
                print("\nExiting help system. Goodbye!\n")
                break
            
            choice_num = int(choice)
            if 1 <= choice_num <= len(topics):
                topic = topics[choice_num - 1]
                print("\n" + "="*70)
                print(help_content[topic])
                print("="*70)
                input("\nPress Enter to continue...")
                print("\n" + "-"*70)
                print("Select another topic or 0 to exit:")
                for i, topic in enumerate(topics, 1):
                    print(f"  {i}. {topic}")
                print("  0. Exit")
                print("-"*70)
            else:
                print("Invalid choice. Please enter a number between 0 and {}.".format(len(topics)))
        
        except ValueError:
            print("Invalid input. Please enter a number.")
        except KeyboardInterrupt:
            print("\n\nExiting help system. Goodbye!\n")
            break


def main():
    """Main function for standalone execution."""
    import sys
    
    if len(sys.argv) > 1:
        # Command-line argument provided
        arg = sys.argv[1].lower()
        
        if arg in ['--help', '-h']:
            print(__doc__)
            return
        
        # Check if requesting specific tab
        if arg.startswith('--tab='):
            tab_name = arg.split('=')[1]
        elif arg == '--tab' and len(sys.argv) > 2:
            tab_name = sys.argv[2]
        else:
            tab_name = arg
        
        # Try to find matching tab
        help_content = get_help_content()
        matched = None
        
        for key in help_content.keys():
            if tab_name.lower() in key.lower().replace(' ', ''):
                matched = key
                break
        
        if matched:
            print("\n" + "="*70)
            print(help_content[matched])
            print("="*70 + "\n")
        else:
            print(f"\nError: Unknown tab '{tab_name}'")
            print("\nAvailable tabs:")
            for key in help_content.keys():
                print(f"  - {key}")
            print("\nUsage: python help.py --tab=<name>")
            print("       python help.py <name>\n")
    else:
        # No arguments - show interactive menu
        display_help_menu()


if __name__ == "__main__":
    main()
