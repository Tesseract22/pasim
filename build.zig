const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const raylib = b.dependency("raylib", .{ .target = target, .optimize = optimize });

    const pasim = b.addModule("pasim", .{
        .optimize = optimize,
        .target = target,
        .root_source_file = b.path("src/main.zig"),
    });

    // raylib 
    pasim.addIncludePath(raylib.path("src"));
    pasim.linkLibrary(raylib.artifact("raylib"));

    // raygui
    pasim.addCSourceFile(.{
        .file = b.path("src/thirdparty/raygui.h"),
        .flags = &.{ "-DRAYGUI_IMPLEMENTATION" },
        .language = .c,
    });
    pasim.addIncludePath(b.path("src/thirdparty"));

    const exe = b.addExecutable(.{
        .name = "pasim",
        .root_module = pasim,
    });

    b.installArtifact(exe);
}
