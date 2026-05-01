
import toga, os, asyncio
from toga.style import Pack
from toga.style.pack import COLUMN, ROW, LEFT,CENTER

class FramesCounter:
    
    async def open_file_dialog(self, widget):
        try:
            self.trajec = await self.main_window.dialog(
                toga.OpenFileDialog("Open file")
            )

            if self.trajec is not None:
                self.textInput_file.value = f"{self.trajec}"
            else:
                await self.main_window.dialog(
                    toga.InfoDialog("Warning", "No file was selected!")
                )
                return

        except ValueError:
            await self.main_window.dialog(
                toga.InfoDialog("Error", "Open file was canceled!")
            )
            return

        self.output_dir = os.path.dirname(self.trajec)

        try:
            with open(self.trajec, "r") as f:
                first_line = f.readline()
                if not first_line:
                    await self.main_window.dialog(
                        toga.InfoDialog("Error", "The TRAJEC.xyz file is empty!")
                    )
                    return
                self.num_atoms = int(first_line.strip().split()[0])
        except Exception as e:
            await self.main_window.dialog(
                toga.InfoDialog("Error", f"Failed to read the file: {e}")
            )
            return

    async def open_geometry_xyz(self,widget):
        self.coords = []
        await self.open_file_dialog(widget)
        with open(self.trajec) as f:
            lines = f.readlines()

        for i in range(2, self.num_atoms + 2):
            atom_data = lines[i].split()
            atom = atom_data[0]
            x, y, z = map(float, atom_data[1:4])
            self.coords.append([atom, x, y, z])

        atom_list = [coord[0] for coord in self.coords]
        self.element_count = {}

        for element in atom_list:
            self.element_count[element] = self.element_count.get(element, 0) + 1
        
        
        update_text = (f"Number of atoms: {self.num_atoms}\n" f"{self.element_count}")
        self.multi_line_text.value = update_text


    async def frames_counter(self, widget):
        """
        Counts the number of frames in TRAJEC.xyz file
        """
        await self.open_file_dialog(widget)
        if not hasattr(self, "trajec") or not self.trajec:
            return

        self.multi_line_text.value = (
            f"\nReading TRAJEC.xyz file ........\n"
        )

        
        frame_count = 0
        total_lines = 0

        try:
            with open(self.trajec, "r") as f:
                while True:
                    title_line = f.readline()
                    if not title_line:
                        break
                    comment_line = f.readline()
                    if not comment_line:
                        break
                    total_lines += 2

                    for _ in range(self.num_atoms):
                        atom_line = f.readline().split() # Read atom line
                        total_lines += 1
                        if not atom_line or len(atom_line) != 4:
                            raise ValueError(f"Error reading atom line {total_lines} in frame {frame_count}")
                            

                    frame_count += 1

                    if (frame_count % 400) == 0:
                        self.progress_label.text = f"{frame_count}"
                        await asyncio.sleep(0)

                self.progress_label.text = f"{frame_count}"
                await asyncio.sleep(0)

            number_of_atoms = f"Number of atoms  -->  {self.num_atoms}\n"
            number_of_lines = f"Number of lines  -->  {total_lines}\n"
            number_of_frames = f"Number of frames -->  {frame_count}\n"
            first_frame = f"  TRAJEC.xyz FIRST FRAME\n"

            with open(self.trajec, "r") as f:
                frame_lines = [f.readline() for _ in range(self.num_atoms + 2)]

            update_text = (
                number_of_atoms
                + number_of_lines
                + number_of_frames
                + first_frame
                + "".join(frame_lines)
            )
            self.multi_line_text.value = update_text

            # self.frame_count = frame_count  # Store for later use
            self.total_frame_number = frame_count  # Store for later use
        except Exception as e:
            await self.main_window.dialog(
                toga.InfoDialog("Error", f"Error reading TRAJEC.xyz file: {e}")
            )
