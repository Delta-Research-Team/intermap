"""
Configuration settings for the InterMap Visualizations app.

"""

# Color scheme for different interaction types
all_interactions_colors = {
    'HBDonor': '#E41A1C',
    'HBAcceptor': '#377EB8',
    'Hydrophobic': '#4DAF4A',
    'VdWContact': '#FF7F00',
    'CloseContacts': '#984EA3',
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

# Error messages
ERROR_MESSAGES = {
    'no_file': 'Please upload a CSV file.',
    'invalid_file': 'Invalid file format. Please upload a CSV file.',
    'no_data': 'No data available for plotting.',
    'processing_error': 'Error processing file: {}',
    'no_interactions': 'No interactions found with the current filters.',
    'plot_error': 'Error generating plot: {}'
}

# CSS styles for UI components
CSS_STYLES = {
    'error_message': """
            .error-message {
                color: #ff5555ff !important;
                font-family: 'Ubuntu Mono', monospace !important;
                padding: 10px !important;
                margin: 10px 0 !important;
                border-radius: 4px !important;
                background-color: rgba(255, 85, 85, 0.1) !important;
                display: block !important;
            }
        """,

    'font_import': """
    @import url('https://fonts.googleapis.com/css2?family=Ubuntu+Mono&display=swap');
    """,

    'custom_controls': """
        /* Slider styles */
        .irs--shiny .irs-bar {
            background: #ff5555ff !important;
            border-top: 1px solid #ff5555ff !important;
            border-bottom: 1px solid #ff5555ff !important;
        }
        
        .irs--shiny .irs-single {
            background: #ff5555ff !important;
            border-color: #ff5555ff !important;
        }

        /* Slider handle */
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

        /* Checkbox styles  */
        .checkbox-group input[type="checkbox"] {
            -webkit-appearance: none !important;
            -moz-appearance: none !important;
            appearance: none !important;
            width: 16px !important;
            height: 16px !important;
            border: 2px solid #ff5555ff !important;
            border-radius: 3px !important;
            background-color: #FEFBF6 !important;
            cursor: pointer !important;
            margin-right: 8px !important;
            position: relative !important;
        }

        .checkbox-group input[type="checkbox"]:checked {
            background-color: #ff5555ff !important;
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
        
        /* Switch color */
        .form-switch .form-check-input:checked {
            background-color: #ff5555ff !important;
            border-color: #ff5555ff !important;
        }
        
        /* Button color */
        .search-button {
            background-color: #ff5555ff !important;
            color: white !important;
            border: none !important;
        }
        
        /* Button hover */
        .search-button:hover {
            background-color: #ff3333ff !important;
        }
    """,

    'welcome_section': """
        .welcome-section {
            display: flex;
            flex-direction: row-reverse;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 20px;
            padding: 20px;
            background: #F7F7F7;
            width: 100%;
            font-family: 'Ubuntu Mono', monospace;
            color: #383747ff;
        }
    """,

    'welcome_text': """
        .welcome-text {
            text-align: left;
            padding-right: 20px;
            font-family: 'Ubuntu Mono', monospace;
            color: #383747ff;
        }
    """,

    'welcome_title': """
        .welcome-title {
            font-size: 36px;
            font-weight: bold;
            margin: 0;
            color: #383747ff;
            font-family: 'Ubuntu Mono', monospace;
        }
    """,

    'welcome_subtitle': """
        .welcome-subtitle {
            font-size: 18px;
            color: #383747ff;
            margin-top: 5px;
            font-family: 'Ubuntu Mono', monospace;
        }
    """,

    'welcome_image': """
        .welcome-image {
            max-width: 250px;
            height: auto;
            margin-left: 1px;
        }
    """,

    'interaction_filter': """
        .interaction-filter {
            background: #f8f9fa;
            padding: 15px;
            border-radius: 8px;
            margin-bottom: 15px;
            font-family: 'Ubuntu Mono', monospace;
            color: #383747ff;
        }
    """,

    'checkbox_group': """
        .checkbox-group {
            max-height: 300px;
            overflow-y: auto;
            padding: 10px;
            border: 1px solid #ddd;
            border-radius: 4px;
            font-family: 'Ubuntu Mono', monospace;
        }
    """,

    'interaction_count': """
        .interaction-count {
            color: #383747ff;
            font-size: 0.9em;
            margin-left: 4px;
            font-family: 'Ubuntu Mono', monospace;
        }
    """,

    'search_container': """
        .search-container {
            display: flex;
            align-items: center;
            gap: 10px;
            margin-bottom: 10px;
            font-family: 'Ubuntu Mono', monospace;
        }
    """,

    'search_button': """
        .search-button {
            background-color: #ff5555ff;
            color: white;
            border: none;
            padding: 6px 12px;
            border-radius: 4px;
            cursor: pointer;
        }
    """,

    'search_button_hover': """
        .search-button:hover {
            background-color: #ff3333ff;
        }
    """,

    'not_found': """
        .not-found {
            color: red;
            margin-top: 5px;
            font-size: 0.9em;
            font-family: 'Ubuntu Mono', monospace;
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
            font-family: 'Ubuntu Mono', monospace;
        }
    """,

    'centered_plot': """
        .centered-plot {
            display: flex;
            flex-direction: column;
            align-items: center;
            width: 100%;
            max-width: 100%;
            font-family: 'Ubuntu Mono', monospace;
        }
    """,

    'plot_output_container': """
        .plot-output-container {
            width: 100%;
            display: flex;
            justify-content: center;
            align-items: center;
            font-family: 'Ubuntu Mono', monospace;
        }
    """
}
