#include "OpenGLUtils.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <tuple>

// Reads all lines of text out of a file, and returns a string
// containing all those lines, separated by newline characters.
std::string readFile(const std::string& filePath)
{
    std::ifstream inStream(filePath);
    std::stringstream buffer;
    buffer << inStream.rdbuf();
    return buffer.str();
}

// Checks whether a shader compiled successfully.
// Returns (true, [info log string with warnings]) if compilation was successful.
// Otherwise, (false, [info log string with warnings/errors]) is returned.
std::pair<bool, std::string> checkShaderCompilation(GLuint shaderNumber)
{
    GLint status;
    glGetShaderiv(shaderNumber, GL_COMPILE_STATUS, &status);

    GLint logSize;
    glGetShaderiv(shaderNumber, GL_INFO_LOG_LENGTH, &logSize);

    std::vector<GLchar> logBuffer(logSize);
    glGetShaderInfoLog(shaderNumber, logSize, nullptr, logBuffer.data());

    if (status == GL_TRUE)
        return { true, logSize == 0 ? "" : logBuffer.data() };
    else
        return { false, logBuffer.data() };
}

// Creates a shader of the given type from the source code in the file at the given path.
GLuint createShaderFromSource(const std::string& sourceFilePath, GLenum shaderType)
{
    const std::string source = readFile(sourceFilePath);
    const char* sourceCharPtr = source.c_str();

    GLuint shader = glCreateShader(shaderType);
    glShaderSource(shader, 1, &sourceCharPtr, nullptr);
    glCompileShader(shader);

    return shader;
}

// Creates a shader program given a file path to a vertex shader source code file,
// a fragment shader source code file, and a std::vector of
// (fragment shader color number, fragment shader out variable) pairs,
// which are passed to glBindFragDataLocation().
std::tuple<GLuint, GLuint, GLuint> createShaderProgram(
    const std::string& vertexSourceFilePath,
    const std::string& fragmentSourceFilePath,
    const std::vector<std::tuple<GLuint, const char*>>& fragDataLocations)
{
    // Create and compile the vertex shader
    GLuint vertexShader = createShaderFromSource(vertexSourceFilePath, GL_VERTEX_SHADER);
    auto shaderInfo = checkShaderCompilation(vertexShader);
    if (!shaderInfo.first)
        printf("Vertex shader failed to compile!\n%s\n", shaderInfo.second.c_str());
    else if (!shaderInfo.second.empty())
        printf("Vertex shader compiled with warnings.\n%s\n", shaderInfo.second.c_str());

    // Create and compile the fragment shader
    GLuint fragmentShader = createShaderFromSource(fragmentSourceFilePath, GL_FRAGMENT_SHADER);
    shaderInfo = checkShaderCompilation(fragmentShader);
    if (!shaderInfo.first)
        printf("Vertex shader failed to compile!\n%s\n", shaderInfo.second.c_str());
    else if (!shaderInfo.second.empty())
        printf("Vertex shader compiled with warnings.\n%s\n", shaderInfo.second.c_str());

    // Link the vertex and fragment shader into a shader program
    GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    for (const auto& pair : fragDataLocations)
        glBindFragDataLocation(shaderProgram, std::get<0>(pair), std::get<1>(pair));
    glLinkProgram(shaderProgram);
    return std::make_tuple(vertexShader, fragmentShader, shaderProgram);
}
