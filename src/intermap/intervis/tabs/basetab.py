# Created by rglez at 7/27/25
"""
BaseTab class for the InterVis application.

This class serves as a base for all tabs in the InterMap application,
providing common functionality and properties.
"""
import pandas as pd


class BaseTab:
    """
    Base class for all tabs in the InterVis application.
    """

    def __init__(self, df, width=600, height=800, **kwargs):
        """
        Initializes the BaseTab with the provided DataFrame and plot dimensions.

        Args:
            df (pd.DataFrame): DataFrame containing interaction data.
            width (int): Width of the plot.
            height (int): Height of the plot.
            **kwargs: Additional keyword arguments for customization.
        """
        self.df = df
        self.width = width
        self.height = height
        self.processed_data = self.process_data()
        self.show_prevalence = kwargs.get('show_prevalence', None)

    def process_data(self):
        """
        Process the DataFrame to prepare it for visualization.

        This method should be overridden in subclasses to implement specific
        data processing logic.
        """
        raise NotImplementedError("Subclasses must implement this method.")

# =============================================================================
#
# =============================================================================

#pickle = '/home/rglez/RoyHub/intermap/data/tutorial-mayank/outputs/ligs-channel_InterMap.pickle'
#cfg = '/home/rglez/RoyHub/intermap/data/tutorial-mayank/outputs/ligs-channel_InterMap.cfg'
#csv_obj = CSVFilter(pickle, cfg)
#master_df = csv_obj.master

#self = HeatMap(master_df, width=800, height=600, show_prevalence=True)
