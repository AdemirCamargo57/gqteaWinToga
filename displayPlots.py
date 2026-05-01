import toga, tempfile
from toga.style import Pack
import matplotlib
matplotlib.use('Agg') # use non-interactive backend to avoid GUI issues
import matplotlib.pyplot as plt

class DisplayPlots():

    font_style = {'color':  'darkred','weight': 'normal','size': 14}
    saved_plot_files = []

    def save_plots(self, k, x ,y, plot_xlabel, plot_ylabel, plot_title):
        temp_filename = tempfile.NamedTemporaryFile(delete=False, suffix=".png", dir=self.output_dir).name
        plt.figure(k)
        plt.plot(x, y,antialiased=True)
        plt.xlabel(plot_xlabel, fontdict = self.font_style)
        plt.ylabel(plot_ylabel, fontdict = self.font_style)
        plt.title(plot_title, fontdict = self.font_style)
        plt.savefig(temp_filename)
        plt.close()
        self.saved_plot_files.append(temp_filename)

    def display_plots(self):

        for plot_filename in self.saved_plot_files:

            # Load the image using Toga's Image class
            plot_image = toga.Image(plot_filename)

            # Create a new window for the plot
            plot_window = toga.Window(
                title=" ",
                size=(700, 600),
            )

            # Create a box to hold the image
            plot_box = toga.Box(style=Pack(flex=1))
            plot_window.content = plot_box

            # Create an ImageView to display the image
            plot_imageview = toga.ImageView(plot_image, style=Pack(flex=1))
            plot_box.add(plot_imageview)

            # Show the window
            plot_window.show()

        # # Clear the list of saved plot files after displaying them
        self.saved_plot_files = []

