"""
Configuration settings for the InterMap Visualizations app.
Contains color schemes, plot settings, and UI styles.

"""

# Color scheme for different interaction types
"""
all_interactions_colors = {
    'HBDonor': '#E41A1C',
    'HBAcceptor': '#377EB8',
    'Hydrophobic': '#4DAF4A',
    'VdWContact': '#FF7F00',
    'CloseContact': '#984EA3',
    'PiStacking': '#F781BF',
    'PiCation': '#A65628',
    'CationPi': '#999999',
    'Anionic': '#B2182B',
    'Cationic': '#2166AC',
    'MetalDonor': '#A8A8A8',
    'MetalAcceptor': '#8C510A',
    'FaceToFace': '#762A83',
    'EdgeToFace': '#1B7837',
    'XBAcceptor': '#5E3C99',
    'XBDonor': '#2B5C8C'
}
"""

all_interactions_colors = {
    'HBDonor': '#FF6B6B',
    'HBAcceptor': '#4ECDC4',
    'Cationic': '#2B95FF',
    'Anionic': '#FF8F40',
    'WaterBridge': '#87CEEB',
    'PiStacking': '#9B6B9E',
    'PiCation': '#BE79BE',
    'PiAnion': '#418641',
    'CationPi': '#D4A5D4',
    'FaceToFace': '#845B87',
    'EdgeToFace': '#AA8BAD',
    'MetalDonor': '#B8B8B8',
    'MetalAcceptor': '#D4AF37',
    'VdWContact': '#D3D3D3',
    'CloseContact': '#E6E6E6',
    'Hydrophobic': '#90A4AE',
    'XBAcceptor': '#66CDAA',
    'XBDonor': '#20B2AA'
}

# Error messages dictionary
ERROR_MESSAGES = {
    'nothing_selected': "Change your filters, your selection is empty",
    'no_file': "Please upload a pickle file.",
    'invalid_file': "Invalid file format. Please upload a pickle file.",
    'no_data': "No data available for plotting.",
    'processing_error': "Error processing file: {}",
    'no_interactions': "No interactions found with current filters.",
    'plot_error': "Error generating plot: {}",
    'no_topology': "No configuration file found in directory.",
    'invalid_sele': "Invalid MDAnalysis selection: {}"
}

# CSS styles for UI components
CSS_STYLES = {
    'file_input_container': """
    .file-input-container {
        margin-bottom: 15px;
    }
    """,

    'file_browse_container': """
    .file-browse-container {
        display: flex;
        align-items: center;
        gap: 10px;
        margin-bottom: 10px;
    }
    """,

    'topology_indicator': """
    .topology-indicator {
        display: inline-flex;
        align-items: center;
        margin-left: 10px;
        cursor: help;
        position: relative;
    }
    
    .topology-indicator i {
        font-size: 14px;
        color: #808080;
        transition: color 0.3s ease;
    }
    
    .topology-indicator.active i {
        color: #4CAF50;
    }
    
    .topology-indicator::after {
        content: attr(title);
        position: absolute;
        bottom: 100%;
        left: 50%;
        transform: translateX(-50%);
        padding: 4px 8px;
        background-color: rgba(0, 0, 0, 0.8);
        color: white;
        border-radius: 4px;
        font-size: 12px;
        white-space: nowrap;
        opacity: 0;
        visibility: hidden;
        transition: opacity 0.3s ease, visibility 0.3s ease;
    }
    
    .topology-indicator:hover::after {
        opacity: 1;
        visibility: visible;
    }
    """,

    'mda_selection_container': """
    .mda-selection-container {
        margin-top: 10px;
        margin-bottom: 15px;
    }
    
    .mda-selection-container input {
        font-family: 'Roboto', sans-serif;
        padding: 8px;
        border: 1px solid #ddd;
        border-radius: 4px;
        background-color: #fff;
        width: 100%;
    }
    
    .mda-selection-container input:focus {
        border-color: #4051b5ff;
        outline: none;
        box-shadow: 0 0 0 2px rgba(255,85,85,0.2);
    }
    
    .mda-selection-container label {
        display: block;
        margin-bottom: 5px;
        font-family: 'Roboto', sans-serif;
        color: #383747ff;
    }
    """,

    'file_upload': """
    .progress-bar {
        background-color: #4051b5ff !important;
    }

    .btn-file {
        background-color: #383747ff !important;
        color: white !important;
        border: none !important;
        padding: 6px 12px !important;
        border-radius: 4px !important;
        cursor: pointer !important;
        transition: background-color 0.3s ease !important;
    }

    .btn-file:hover {
        background-color: #4051b5ff !important;
    }

    .progress {
        border-radius: 4px !important;
        margin-top: 5px !important;
        background-color: #f0f0f0 !important;
    }
    """,

    'error_message': """
    .error-message {
        color: #4051b5ff !important;
        font-family: 'Roboto', sans-serif !important;
        padding: 10px !important;
        margin: 10px 0 !important;
        border-radius: 4px !important;
        background-color: rgba(255, 85, 85, 0.1) !important;
        display: block !important;
        border-left: 4px solid #4051b5ff !important;
    }
    """,

    'font_import': """
    @import url('https://fonts.googleapis.com/css2?family=Roboto:wght@300;400;500;700&display=swap');
    """,

    'custom_controls': """
    .irs--shiny .irs-bar {
        background: #4051b5ff !important;
        border-top: 1px solid #4051b5ff !important;
        border-bottom: 1px solid #4051b5ff !important;
    }

    .irs--shiny .irs-single {
        background: #4051b5ff !important;
        border-color: #4051b5ff !important;
    }

    .irs--shiny .irs-handle {
        border: 3px solid #383747ff !important;
        background-color: #383747ff !important;
    }

    .irs--shiny .irs-handle > i:first-child {
        background-color: #383747ff !important;
    }

    .irs--shiny .irs-handle.state_hover > i:first-child,
    .irs--shiny .irs-handle:hover > i:first-child {
        background-color: #383747ff !important;
    }

    .checkbox-group input[type="checkbox"] {
        -webkit-appearance: none !important;
        -moz-appearance: none !important;
        appearance: none !important;
        width: 16px !important;
        height: 16px !important;
        border: 2px solid #4051b5ff !important;
        border-radius: 3px !important;
        background-color: #FEFBF6 !important;
        cursor: pointer !important;
        margin-right: 8px !important;
        position: relative !important;
        transition: all 0.3s ease !important;
    }

    .checkbox-group input[type="checkbox"]:checked {
        background-color: #4051b5ff !important;
    }

    .checkbox-group input[type="checkbox"]:checked::after {
        content: '✓' !important;
        color: #FEFBF6 !important;
        position: absolute !important;
        left: 50% !important;
        top: 50% !important;
        transform: translate(-50%, -50%) !important;
        font-size: 12px !important;
    }

    .form-switch .form-check-input:checked {
        background-color: #4051b5ff !important;
        border-color: #4051b5ff !important;
    }

    .search-button {
        background-color: #4051b5ff !important;
        color: white !important;
        border: none !important;
        padding: 6px 12px !important;
        border-radius: 4px !important;
        cursor: pointer !important;
        transition: background-color 0.3s ease !important;
        font-family: 'Roboto', sans-serif !important;
    }

    .search-button:hover {
        background-color: #ff3333ff !important;
    }

    .form-control:focus {
        border-color: #4051b5ff !important;
        box-shadow: 0 0 0 0.2rem #4051b5ff !important;
    }
    
    .welcome-section {
        display: flex;
        flex-direction: row-reverse;
        justify-content: space-between;
        align-items: center;
        margin: 0;
        padding: 5px;
        background: #4051b5ff;
        width: 100vw;
        font-family: 'Roboto', sans-serif;
        color: #383747ff;
        position: fixed;          
        top: 0;                  
        left: 0;                 
        right: 0;                
        margin-left: 0vw;
        margin-right: 0vw;
        z-index: 999;
        box-shadow: 0 3px 6px rgba(0, 0, 0, 0.16), 0 3px 6px rgba(0, 0, 0, 0.23);
    }
    
        .row {
            margin-top: 120px; 
        }
    """


    ,

    'welcome_text': """
    .welcome-text {
        text-align: right;
        padding-right: 20px;
        font-family: 'Roboto', sans-serif;
        color: white;
    }
    """,

    'welcome_title': """
    .welcome-title {
        font-size: 36px;
        font-weight: bold;
        margin: 0;
        color: white;
        font-family: 'Roboto', sans-serif;  
    }
    """,

    'welcome_subtitle': """
    .welcome-subtitle {
        font-size: 20px;
        color: white;
        margin-top: 5px;
        font-family: 'Roboto', sans-serif;
        opacity: 0.8;
    }
    """,

    'welcome_image': """
    .welcome-image {
        max-width: 250px;
        height: auto;
        margin-left: 15px;
        border-radius: 8px;
        box-shadow: none;
    }
    """,

    'interaction_filter': """
    .interaction-filter {
        background: #f8f9fa;
        padding: 15px;
        border-radius: 8px;
        margin-bottom: 15px;
        font-family: 'Roboto', sans-serif;
        color: #383747ff;
        box-shadow: 0 2px 4px rgba(0,0,0,0.05);
    }
    """,

    'checkbox_group': """
    .checkbox-group {
        max-height: 300px;
        overflow-y: auto;
        padding: 10px;
        border: 1px solid #ddd;
        border-radius: 4px;
        font-family: 'Roboto', sans-serif;
        scrollbar-width: thin;
        scrollbar-color: #4051b5ff #f0f0f0;
    }

    .checkbox-group::-webkit-scrollbar {
        width: 8px;
    }

    .checkbox-group::-webkit-scrollbar-track {
        background: #f0f0f0;
        border-radius: 4px;
    }

    .checkbox-group::-webkit-scrollbar-thumb {
        background-color: #4051b5ff;
        border-radius: 4px;
    }
    """,

    'interaction_count': """
    .interaction-count {
        color: #383747ff;
        font-size: 0.9em;
        margin-left: 4px;
        font-family: 'Roboto', sans-serif;
        opacity: 0.8;
    }
    """,

    'select_all_checkbox': """
    .interaction-filter .select-all-container {
        padding: 4px !important;
        margin-bottom: 5px !important;
        background-color: transparent !important;
        border-radius: 4px !important;
    }

    .interaction-filter .select-all-container input[type="checkbox"] {
        -webkit-appearance: none !important;
        -moz-appearance: none !important;
        appearance: none !important;
        width: 16px !important;
        height: 16px !important;
        border: 2px solid #ff5555ff !important; 
        border-radius: 50% !important;
        background-color: white !important;
        cursor: pointer !important;
        margin-right: 8px !important;
        position: relative !important;
        transition: all 0.3s ease !important;
    }

    .interaction-filter .select-all-container input[type="checkbox"]:checked {
        background-color: #ff5555ff !important;  
    }

    .interaction-filter .select-all-container input[type="checkbox"]:checked::after {
        content: '✓' !important;
        color: white !important;
        position: absolute !important;
        left: 50% !important;
        top: 50% !important;
        transform: translate(-50%, -50%) !important;
        font-size: 12px !important;
    }

    .interaction-filter .select-all-container label {
        font-family: 'Roboto', sans-serif !important;
        font-size: 0.9em !important;
        color: black !important;  
    }
    """,

    'plot_container': """
    .plot-container {
        display: flex;
        justify-content: center;
        align-items: center;
        width: 100%;
        margin: 20px auto;
        padding: 0 15px;
        font-family: 'Roboto', sans-serif;
    }
    """,

    'centered_plot': """
    .centered-plot {
        display: flex;
        flex-direction: column;
        align-items: center;
        width: 100%;
        max-width: 100%;
        font-family: 'Roboto', sans-serif;
    }
    """,

    'plot_output_container': """
    .plot-output-container {
        width: 100%;
        display: flex;
        justify-content: center;
        align-items: center;
        font-family: 'Roboto', sans-serif;
        background: white;
        border-radius: 8px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.05);
        padding: 15px;
        margin-bottom: 20px;
    }
    """,

    'tabs': """
        .nav-tabs {
            border-bottom: none;
            display: flex;
            justify-content: center;
            gap: 10px;
            padding: 10px 0;
            background-color: #ffffff;
            margin-bottom: 20px;
        }

        .nav-tabs .nav-link {
            position: relative;
            background-color: #f0f0f0;
            border: none;
            border-radius: 6px;
            color: #333333;
            padding: 12px 24px;
            margin: 0;
            transition: all 0.3s ease;
            font-family: 'Roboto', sans-serif;
            font-size: 14px;
            font-weight: 500;
            text-align: center;
            min-width: 150px;
        }

        .nav-tabs .nav-link:hover {
            background-color: #4051b5ff;  
            color: white;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }

        .nav-tabs .nav-link.active {
            background-color: #4051b5ff !important;  /* Cambiado a azul */
            color: white !important;  /* Cambiado a blanco */
            font-weight: bold;
            border: none;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }

        .nav-tabs .nav-link.active:hover {
            background-color: #4051b5ff;  
            color: white;
        }

        .custom-tabs-container {
            width: 100%;
            margin-bottom: 20px;
        }
        
        .control-group {
            background: white;
            padding: 12px;
            border-radius: 6px;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        }
        
        .control-group h5 {
            font-size: 0.9em;
            font-weight: 500;
            margin-bottom: 8px;
        }
        
        .form-group {
            margin-bottom: 10px;
        }
        
        .form-group label {
            font-size: 0.85em;
            color: #666;
        }
        
        .slider-container {
            padding: 0 10px;
        }
        
        .custom-checkbox {
            margin-bottom: 8px;
        }
        """

}
