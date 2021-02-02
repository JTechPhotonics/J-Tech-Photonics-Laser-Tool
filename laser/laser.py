import os

from lxml import etree
from inkex import EffectExtension, Boolean

from svg_to_gcode.svg_parser import parse_string
from svg_to_gcode.geometry import LineSegmentChain
from svg_to_gcode.compiler import Compiler, interfaces

inx_filename = "laser.inx"


class GcodeExtension(EffectExtension):

    def __init__(self):
        EffectExtension.__init__(self)

    def effect(self):
        interface = {"gcode": interfaces.Gcode, "gcode_fan": interfaces.FanControlledGcode}[self.options.interface]
        movement_speed = self.options.travel_speed
        cutting_speed = self.options.laser_speed
        pass_depth = self.options.pass_depth
        passes = self.options.passes
        output_path = os.path.join(self.options.directory, self.options.filename)
        filename_suffix = self.options.filename_suffix

        if filename_suffix:
            filename, extension = output_path.split('.')

            n = 1
            while os.path.isfile(output_path):
                output_path = filename + str(n) + '.' + extension
                n += 1

        gcode_compiler = Compiler(interface, movement_speed, cutting_speed, pass_depth)

        curves = parse_string(etree.tostring(self.document))  # Parse an svg file into geometric curves

        gcode_compiler.append_curves(curves)
        gcode_compiler.compile_to_file(output_path, passes=passes)

        return self.document

    def draw_debug_lines(self, svg, curves):
        name_space = 'http://www.w3.org/2000/svg'

        tolerance = self.options.approximation_tolerance
        unit = self.options.unit

        tree = etree.ElementTree()
        tree.fromstring(svg)

        root = tree.getroot()

        group = etree.Element("{%s}g" % name_space)
        group.set("id", "debug_layer")
        group.set("inkscape:groupmode", "layer")
        group.set("inkscape:label", "debug laser output layer")

        for curve in curves:
            approximation = LineSegmentChain.line_segment_approximation(curve)

            path = etree.Element("{%s}path" % name_space)
            path.set("d", approximation.to_svg_path(wrapped=False))
            path.set("fill", "none")
            path.set("stroke", "red")
            path.set("stroke-width", f"{tolerance / 2}{unit}")

            group.append(path)

        root.append(group)

        return tree

    def add_arguments(self, arg_parser):
        arguments = self.read_arguments()

        for arg in arguments:
            arg_parser.add_argument("--" + arg["name"], type=arg["type"], dest=arg["name"])

    @staticmethod
    def read_arguments():
        root = etree.parse(inx_filename).getroot()

        arguments = []  # [{name, type, ...}]
        namespace = "http://www.inkscape.org/namespace/inkscape/extension"
        for arg in root.iter("{%s}param" % namespace):

            name = arg.attrib["name"]

            arg_type = arg.attrib["type"]

            if arg_type in ["description"]:
                continue

            types = {"int": int, "float": float, "boolean": Boolean, "string": str, "enum": str}

            arguments.append({"name": name, "type": types[arg_type]})

        return arguments


if __name__ == '__main__':
    effect = GcodeExtension()
    effect.run()
