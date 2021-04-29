import os

from lxml import etree
from xml.etree import ElementTree as xml_tree
from inkex import EffectExtension, Boolean

from svg_to_gcode.svg_parser import parse_root, Transformation, debug_methods
from svg_to_gcode.geometry import LineSegmentChain
from svg_to_gcode.compiler import Compiler, interfaces
from svg_to_gcode.formulas import linear_map
from svg_to_gcode import TOLERANCES

svg_name_space = "http://www.w3.org/2000/svg"
inkscape_name_space = "http://www.inkscape.org/namespaces/inkscape"
sodipodi_name_space = "http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"

inx_filename = "laser.inx"


def generate_custom_interface(laser_off_command, laser_power_command, laser_power_range):
    """Wrapper function for generating a Gcode interface with a custom laser power command"""

    class CustomInterface(interfaces.Gcode):
        """A Gcode interface with a custom laser power command"""

        def __init__(self):
            super().__init__()

        def laser_off(self):
            return f"{laser_off_command};"

        def set_laser_power(self, power):
            if power < 0 or power > 1:
                raise ValueError(f"{power} is out of bounds. Laser power must be given between 0 and 1. "
                                 f"The interface will scale it correctly.")

            return f"{laser_power_command} S{linear_map(0, int(laser_power_range), power)};"

    return CustomInterface


class GcodeExtension(EffectExtension):
    """Inkscape Effect Extension."""

    def __init__(self):
        EffectExtension.__init__(self)

    def effect(self):
        """Takes the SVG from Inkscape, generates gcode, returns the SVG after adding debug lines."""

        root = self.document.getroot()

        approximation_tolerance = float(self.options.approximation_tolerance.replace(',', '.'))

        output_path = os.path.join(self.options.directory, self.options.filename)
        if self.options.filename_suffix:
            filename, extension = output_path.split('.')

            n = 1
            while os.path.isfile(output_path):
                output_path = filename + str(n) + '.' + extension
                n += 1

        header = None
        if os.path.isfile(self.options.header_path):
            logger.debug(F"going to read{self.options.header_path}")
            with open(self.options.header_path, 'r') as header_file:
                header = header_file.read().splitlines()
                logger.debug(F"This is my header: >>>{header}<<<")
        elif self.options.header_path != os.getcwd():
            self.debug(f"Header file does not exist at {self.options.header_path}")

        if self.options.set_z_axis_start_pos:
            temp = F"G1 Z{self.options.z_axis_start};"
            if header is None:
                header = [temp]
            else:
                header.append(temp)

        footer = None
        if os.path.isfile(self.options.footer_path):
            with open(self.options.footer_path, 'r') as footer_file:
                footer = footer_file.read().splitlines()
        elif self.options.footer_path != os.getcwd():
            self.debug(f"Footer file does not exist at {self.options.footer_path}")

        if self.options.move_to_zero_at_end:
            temp = "G1 X0.0 Y0.0 Z0.0;"
            if footer is None:
                footer = [temp]
            else:
                footer.append(temp)

        # Generate gcode
        self.clear_debug()

        TOLERANCES["approximation"] = approximation_tolerance
        custom_interface = generate_custom_interface(self.options.laser_off_command, self.options.laser_power_command,
                                                     self.options.laser_power_range)

        gcode_compiler = Compiler(custom_interface, self.options.travel_speed, self.options.cutting_speed,
                                  self.options.pass_depth, dwell_time=self.options.dwell_time, custom_header=header,
                                  custom_footer=footer, unit=self.options.unit)

        transformation = Transformation()

        transformation.add_translation(self.options.horizontal_offset, self.options.vertical_offset)
        transformation.add_scale(self.options.scaling_factor)

        if self.options.machine_origin == "center":
            transformation.add_translation(-self.options.bed_width / 2, self.options.bed_height / 2)

        curves = parse_root(root, transform_origin=not self.options.invert_y_axis, root_transformation=transformation,
                            canvas_height=self.options.bed_height)

        gcode_compiler.append_curves(curves)
        gcode_compiler.compile_to_file(output_path, passes=self.options.passes)

        # Generate debug lines
        if self.options.draw_debug:
            self.draw_debug_traces(curves)
            self.draw_unit_reference()
            self.select_non_debug_layer()

        return self.document

    def draw_debug_traces(self, curves):
        """Traces arrows over all parsed paths"""

        root = self.document.getroot()
        origin = self.options.machine_origin
        bed_width = self.options.bed_width
        bed_height = self.options.bed_height

        group = etree.Element("{%s}g" % svg_name_space)
        group.set("id", "debug_traces")
        group.set("{%s}groupmode" % inkscape_name_space, "layer")
        group.set("{%s}label" % inkscape_name_space, "debug traces")

        group.append(
            etree.fromstring(xml_tree.tostring(debug_methods.arrow_defs(arrow_scale=self.options.debug_arrow_scale))))

        for curve in curves:
            approximation = LineSegmentChain.line_segment_approximation(curve)

            change_origin = Transformation()

            if not self.options.invert_y_axis:
                change_origin.add_scale(1, -1)
                change_origin.add_translation(0, -bed_height)

            if origin == "center":
                change_origin.add_translation(bed_width / 2, bed_height / 2)



            path_string = xml_tree.tostring(
                debug_methods.to_svg_path(approximation, color="red", opacity="0.5", stroke_width=f"{self.options.debug_line_width}px",
                                          transformation=change_origin, draw_arrows=True)
            )

            group.append(etree.fromstring(path_string))

        root.append(group)

    def draw_unit_reference(self):
        """Draws reference points to mark the bed's four corners"""
        root = self.document.getroot()
        unit = self.options.unit
        origin = self.options.machine_origin
        bed_width = self.options.bed_width
        bed_height = self.options.bed_height

        group = etree.Element("{%s}g" % svg_name_space)
        group.set("id", "debug_references")
        group.set("{%s}groupmode" % inkscape_name_space, "layer")
        group.set("{%s}label" % inkscape_name_space, "debug reference points")

        reference_points_svg = [(0, 0), (0, bed_height), (bed_width, 0), (bed_width, bed_height)]
        reference_points_gcode = {
            "bottom-left": [(0, bed_height), (0, 0), (bed_width, bed_height), (bed_width, 0)],
            "top-left": [(0, 0), (0, bed_height), (bed_width, 0), (bed_width, bed_height)],
            "center": [(-bed_width / 2, bed_height / 2), (-bed_width / 2, -bed_height / 2),
                       (bed_width / 2, bed_height / 2),
                       (bed_width / 2, -bed_height / 2)]
        }[origin]
        for i, (x, y) in enumerate(reference_points_svg):
            reference_point = etree.Element("{%s}g" % svg_name_space)

            stroke_width = 2
            size = 7

            x_direction = -1 if x > 0 else 1
            plus_sign = etree.Element("{%s}g" % svg_name_space)
            horizontal = etree.Element("{%s}line" % svg_name_space)
            horizontal.set("x1", str(x - x_direction * stroke_width / 2))
            horizontal.set("y1", str(y))
            horizontal.set("x2", str(x + x_direction * size))
            horizontal.set("y2", str(y))
            horizontal.set("style", f"stroke:black;stroke-width:{stroke_width}")
            plus_sign.append(horizontal)

            y_direction = -1 if y > 0 else 1
            vertical = etree.Element("{%s}line" % svg_name_space)
            vertical.set("x1", str(x))
            vertical.set("y1", str(y + stroke_width / 2))
            vertical.set("x2", str(x))
            vertical.set("y2", str(y + y_direction * size))
            vertical.set("style", f"stroke:black;stroke-width:{stroke_width}")
            plus_sign.append(vertical)

            reference_point.append(plus_sign)

            text_box = etree.Element("{%s}text" % svg_name_space)
            text_box.set("x", str(x - 28))
            text_box.set("y", str(y - (y <= 0) * 6 + (y > 0) * 9))
            text_box.set("font-size", "6")
            text_box.text = f"{reference_points_gcode[i][0]}{unit}, {reference_points_gcode[i][1]}{unit}"
            reference_point.append(text_box)

            group.append(reference_point)

        root.append(group)

    def select_non_debug_layer(self):
        """
        Select content_layer and create one if it doesn't exist. This helps stop the user from accidentally placing new
        objects in debug layers.
        """

        root = self.document.getroot()

        unique_id = "layer89324"
        content_layer = root.find("{%s}g[@id='%s']" % (svg_name_space, unique_id))

        if content_layer is None:
            content_layer = etree.Element("{%s}g" % svg_name_space)
            content_layer.set("id", unique_id)
            content_layer.set("{%s}groupmode" % inkscape_name_space, "layer")
            content_layer.set("{%s}label" % inkscape_name_space, "content layer")

        sodipodi = root.find("{%s}namedview" % sodipodi_name_space)
        if sodipodi is not None:
            sodipodi.set("{%s}current-layer" % inkscape_name_space, unique_id)

        root.append(content_layer)

    def clear_debug(self):
        """Removes debug groups. Used before parsing paths for gcode."""
        root = self.document.getroot()

        debug_traces = root.find("{%s}g[@id='debug_traces']" % svg_name_space)
        debug_references = root.find("{%s}g[@id='debug_references']" % svg_name_space)

        if debug_traces is not None:
            root.remove(debug_traces)

        if debug_references is not None:
            root.remove(debug_references)

    def add_arguments(self, arg_parser):
        """Tell inkscape what arguments to stick in self.options (behind the hood it's more complicated, see docs)"""
        arguments = self.read_arguments()

        for arg in arguments:
            arg_parser.add_argument("--" + arg["name"], type=arg["type"], dest=arg["name"])

    @staticmethod
    def read_arguments():
        """
        This method reads arguments off of the inx file so you don't have to explicitly declare them in self.add_arguments()
        """
        root = etree.parse(inx_filename).getroot()

        arguments = []  # [{name, type, ...}]
        namespace = "http://www.inkscape.org/namespace/inkscape/extension"
        for arg in root.iter("{%s}param" % namespace):

            name = arg.attrib["name"]

            arg_type = arg.attrib["type"]

            if arg_type in ["description", "notebook"]:
                continue

            types = {"int": int, "float": float, "boolean": Boolean, "string": str, "enum": str, "path": str}

            arguments.append({"name": name, "type": types[arg_type]})

        if next(root.iter("{%s}page" % namespace)) is not None:
            arguments.append({"name": "tabs", "type": str})

        return arguments


if __name__ == '__main__':
    effect = GcodeExtension()
    effect.run()
