#pragma once

#include <string>
#include <vector>
#include <GL/glew.h>

// Reads all lines of text out of a file, and returns a string
// containing all those lines, separated by newline characters.
std::string readFile(const std::string& filePath);

// Checks whether a shader compiled successfully.
// Returns (true, [info log string with warnings]) if compilation was successful.
// Otherwise, (false, [info log string with warnings/errors]) is returned.
std::pair<bool, std::string> checkShaderCompilation(GLuint shaderNumber);

// Creates a shader of the given type from the source code in the file at the given path.
GLuint createShaderFromSource(const std::string& sourceFilePath, GLenum shaderType);

// Creates a shader program given a file path to a vertex shader source code file,
// a fragment shader source code file, and a std::vector of
// (fragment shader color number, fragment shader out variable) pairs,
// which are passed to glBindFragDataLocation().
// Returns a (vertex shader ID, fragment shader ID, program ID) tuple.
std::tuple<GLuint, GLuint, GLuint> createShaderProgram(
    const std::string& vertexSourceFilePath,
    const std::string& fragmentSourceFilePath,
    const std::vector<std::tuple<GLuint, const char*>>& fragDataLocations);