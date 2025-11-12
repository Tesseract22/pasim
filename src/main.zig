const std = @import("std");

const c = @cImport({
    @cInclude("raylib.h");
    @cInclude("rlgl.h");
});

const w_HEIGHT = 1080;
const w_WIDTH = 1920;
const w_RATIO: comptime_float = 
    @as(comptime_float, @floatFromInt(w_WIDTH)) / 
    @as(comptime_float, @floatFromInt(w_HEIGHT));

const r_VS_PATH = "resources/base.vs";
const r_FS_PATH = "resources/base.fs";
const r_HDR_FS_PATH = "resources/hdr.fs";

const PasimF = f32;

 // normalized world coordinate
 // (0, 1) is the top of screen
const WorldCoord =  @Vector(2, PasimF);
const ScreenCoord =  @Vector(2, PasimF);
const Particle = struct {
    pos: WorldCoord,
    spd: @Vector(2, PasimF),
    mass: f32 = 1,
    kind: u8,
};

var particles: [2000]Particle = undefined;

fn splatv2(f: PasimF) @Vector(2, PasimF) {
    return @splat(f);
}

fn coord2rlcoord(p: WorldCoord) c.Vector2 {
    return .{ .x = p[0], .y = p[1] };
}

fn w2s(p: WorldCoord) ScreenCoord {
    const v = (p + splatv2(1)) / splatv2(2) * splatv2(w_HEIGHT); // this is from 0-1
    return .{ v[0]+(w_WIDTH-w_HEIGHT)/2, w_HEIGHT-v[1] };
}

fn dist2(a: WorldCoord, b: WorldCoord) PasimF {
    const d = a - b;
    return @reduce(.Add, d*d);
}

fn dist(a: WorldCoord, b: WorldCoord) PasimF {
    return @sqrt(dist2(a, b));
}

var particle_default_texture: c.Texture2D = undefined;

const g_PARTICLE_RADIUS: f32 = 0.005;
const g_PARTICLE_TEXTURE_SIZE: f32 = 100;
const g_PARTICLE_TEXTURE_RADIUS: f32 = g_PARTICLE_RADIUS / (g_PARTICLE_TEXTURE_SIZE / w_WIDTH);


// Assume already in shader mode
fn DrawParticle(p: Particle) void {
    const rgb = particle_colors.items[p.kind];
    const sreen_coord = w2s(p.pos);
    // std.log.debug("screen_coord: {}", .{sreen_coord});
    c.DrawTextureEx(
        particle_default_texture,
        .{ .x = sreen_coord[0] - g_PARTICLE_TEXTURE_SIZE/2, .y = sreen_coord[1] - g_PARTICLE_TEXTURE_SIZE/2 },
        0, g_PARTICLE_TEXTURE_SIZE * p.mass,
        rgb);
}

const ForceConfig = struct {
    radius: f32,
    strength: f32, // >0 -> repel, <0 -> attract
};

// const particle_kind = 3;
// const force_config_matrix = [particle_kind][particle_kind]ForceConfig {
//     .{
//         .{ .radius = 0.07, .strength = 0.01 },
//         .{ .radius = 0.09, .strength = -0.04 },
//         .{ .radius = 0.02, .strength = -0.01 },
//     },
//     .{
//         .{ .radius = 0.07, .strength = 0.02 },
//         .{ .radius = 0.10, .strength = -0.03 },
//         .{ .radius = 0.10, .strength = -0.04 },
//     },
//     .{
//         .{ .radius = 0.05, .strength = 0.02 },
//         .{ .radius = 0.05, .strength = -0.03 },
//         .{ .radius = 0.05, .strength = -0.01 },
//     }
//
// };

var particle_kind: u8 = 0;
var particle_force_configs = std.ArrayList(ForceConfig).empty;
var particle_colors = std.ArrayList(c.Color).empty;
const collision_cfg = ForceConfig {
    .radius = 0.02,
    .strength = 7.0,
};
const particle_drag = 5.0;
const random_force_radius_max = 0.1;
const random_force_radius_min = 0.05;
const random_force_strength_max = 0.3;
const random_force_strength_min = 0.01;

fn float_range(random: std.Random, min: f32, max: f32) f32 {
    return random.float(f32) * (max - min) + min;
}

fn random_sign(random: std.Random) f32 {
    return if (random.boolean()) 1 else -1;
}

fn randomize_config(random: std.Random, a: std.mem.Allocator) void {
    particle_force_configs.clearRetainingCapacity();
    particle_colors.clearRetainingCapacity();
    // particle_kind = random.int(u8) % 5 + 2;
    particle_kind = 3;
    for (0..particle_kind) |i| {
        for (0..particle_kind) |j| {
            _ = i;
            _ = j;
            particle_force_configs.append(a, 
                .{
                    .radius = float_range(random, random_force_radius_min, random_force_strength_max),
                    .strength = random_sign(random) * float_range(random, random_force_strength_min, random_force_strength_max),
                })
            catch unreachable;
        }
        const hsv =  c.Vector3 { .x = random.float(f32) * 360, .y = float_range(random, 0.5, 1), .z = float_range(random, 0.5, 1) };
        const rgb = c.ColorFromHSV(hsv.x, hsv.y, hsv.z);
        particle_colors.append(a, rgb) catch unreachable;
    }
}

fn generate_particle(random: std.Random) void {
    // const init_vel_mul = 0.1;
    const init_vel_mul = 0;
    for (&particles) |*p| {
        p.pos = .{ (random.float(f32)-0.5)*2*w_RATIO, (random.float(f32)-0.5)*2 };
        p.spd = .{ (random.float(f32)-0.5)*init_vel_mul, (random.float(f32)-0.5)*init_vel_mul };
        // p.mass = (random.float(f32) + 0.5) * 2;
        p.mass = 1;
        p.kind = random.int(u8) % particle_kind;
    }
}

fn linear_force(cfg: ForceConfig, d: f32) f32 {
    const f = cfg.strength * @max(0, (cfg.radius-@abs(d)) / cfg.radius);
    return f;
}

fn update_game(dt: f32) void {
    for (0..particles.len) |i| {
        const a = &particles[i];
        for (i+1..particles.len) |j| {
            const b = &particles[j];
            const l = a.pos - b.pos;
            const d = dist(a.pos, b.pos);
            const unit_l = if (d == 0) splatv2(0) else l / splatv2(d);
            // if (d > collision_max_dist) continue;
            const collision_force = linear_force(collision_cfg, d) * a.mass * b.mass;

            const interact_cfg_ab = particle_force_configs.items[a.kind * particle_kind + b.kind];
            const interact_cfg_ba = particle_force_configs.items[b.kind * particle_kind + a.kind];
            const interact_force_ab = linear_force(interact_cfg_ab, d) * a.mass * b.mass;
            const interact_force_ba = linear_force(interact_cfg_ba, d) * a.mass * b.mass;
            
            // c.DrawLineEx(coord2rlcoord(w2s(a.pos)), coord2rlcoord(w2s(b.pos)), 1, c.ORANGE);
         
            a.spd += splatv2(dt*interact_force_ab/a.mass) * unit_l;
            b.spd -= splatv2(dt*interact_force_ba/b.mass) * unit_l;

            a.spd += splatv2(dt*collision_force/a.mass) * unit_l;
            b.spd -= splatv2(dt*collision_force/b.mass) * unit_l;

        }
    }

    for (0..particles.len) |i| {
        const a = &particles[i];
        a.pos += a.spd * splatv2(dt);
        if (a.pos[1] < -1) {
            a.pos[1] = -1;
            a.spd[1] = -a.spd[1];
        }
        if (a.pos[0] < -w_RATIO) {
            a.pos[0] = -w_RATIO;
            a.spd[0] = -a.spd[0];
        }
        if (a.pos[1] > 1.0) {
            a.pos[1] = 1.0;
            a.spd[1] = -a.spd[1];
        }
        if (a.pos[0] > w_RATIO) {
            a.pos[0] = w_RATIO;
            a.spd[0] = -a.spd[0];
        }
        a.spd[0] *= @exp(-particle_drag*dt*@abs(a.spd[0]));
        a.spd[1] *= @exp(-particle_drag*dt*@abs(a.spd[1]));
        // a.spd -= a.spd * splatv2(particle_drag * dt / a.mass);
    }
}

pub fn FixedRingBuffer(comptime T: type, comptime size: u32) type {
    return struct {
        const Self = @This();
        data: [size]T = undefined,
        active: std.StaticBitSet(size) = std.StaticBitSet(size).initEmpty(),
        head: u32 = 0,
        tail: u32 = 0,
        count: u32 = 0,
        len: u32 = size,

        pub fn init(len: u32) Self {
            std.debug.assert(len <= size and len != 0);
            return .{ .len = len };
        }

        pub fn init_with_value(len: u32, val: T) Self {
            std.debug.assert(len <= size and len != 0);
            var self = Self { .len = len };
            @memset(&self.data, val);
            return self;

        }

        pub fn push(self: *Self, el: T) void {
            if (self.count == self.len) {
                self.head = (self.head + 1) % self.len;
            } else {
                self.count += 1;
            }
            self.data[self.tail] = el;
            self.active.set(self.tail);
            self.tail = (self.tail + 1) % self.len;
            
        }

        pub fn avg(self: Self) T {
            var sum: T = 0;
            for (0..size) |i| {
                sum += self.data[i];
            }
            return @divTrunc(sum, @as(T, @intCast(size)));
        }

        pub fn remove(self: *Self, i: u32) void {
            self.active.unset(i);
        }

        pub fn at(self: Self, i: u32) T {
            return self.data[(self.head + i) % self.len];
        }

        pub fn last(self: *Self) *T {
            return &self.data[self.tail];
        }

        pub fn is_full(self: Self) bool {
            return self.len == self.count;
        }

        pub fn clear(self: *Self, el: T) void {
            self.head = 0;
            self.tail = 0;
            self.count = 0;
            self.active = std.StaticBitSet(size).initEmpty();
            @memset(&self.data, el);
        }
    };
}


fn LoadRenderTextureHDR(width: c_int, height: c_int) c.RenderTexture2D
{
    var target: c.RenderTexture2D = .{};

    target.id = c.rlLoadFramebuffer(); // Load an empty framebuffer

    if (target.id > 0)
    {
        c.rlEnableFramebuffer(target.id);

        // Create color texture (default to RGBA)
        target.texture.id = c.rlLoadTexture(null, width, height, c.PIXELFORMAT_UNCOMPRESSED_R16G16B16A16, 1);
        target.texture.width = width;
        target.texture.height = height;
        target.texture.format = c.PIXELFORMAT_UNCOMPRESSED_R16G16B16A16;
        target.texture.mipmaps = 1;

        // Create depth renderbuffer/texture
        target.depth.id = c.rlLoadTextureDepth(width, height, true);
        target.depth.width = width;
        target.depth.height = height;
        target.depth.format = 19;       //DEPTH_COMPONENT_24BIT?
        target.depth.mipmaps = 1;

        // Attach color texture and depth renderbuffer/texture to FBO
        c.rlFramebufferAttach(target.id, target.texture.id, c.RL_ATTACHMENT_COLOR_CHANNEL0, c.RL_ATTACHMENT_TEXTURE2D, 0);
        c.rlFramebufferAttach(target.id, target.depth.id, c.RL_ATTACHMENT_DEPTH, c.RL_ATTACHMENT_RENDERBUFFER, 0);

        // Check if fbo is complete with attachments (valid)
        if (c.rlFramebufferComplete(target.id)) std.log.err("FBO: [ID {}] Framebuffer object created successfully", .{target.id});

        c.rlDisableFramebuffer();
    }
    else std.log.err("FBO: Framebuffer object can not be created", .{});

    return target;
}

pub fn main() !void {
    std.log.debug("opengl version: {s}", .{ c.RLGL_VERSION });
    c.InitWindow(w_WIDTH, w_HEIGHT, "Pasim"); 
    particle_default_texture = c.Texture2D {
        .id = c.rlGetTextureIdDefault(),
        .height = 1,
        .width = 1,
        .mipmaps = 1,
        .format = c.PIXELFORMAT_UNCOMPRESSED_R8G8B8A8,
    };

    var rand = std.Random.Xoroshiro128.init(@bitCast(std.time.microTimestamp()));
    const random = rand.random();



    const shader = c.LoadShader(null, r_FS_PATH);
    const shader_radius_uniform = c.GetShaderLocation(shader, "radius");
    c.SetShaderValue(shader, shader_radius_uniform, @ptrCast(&g_PARTICLE_TEXTURE_RADIUS), c.SHADER_UNIFORM_FLOAT);
    const HDRBuffer = LoadRenderTextureHDR(w_WIDTH, w_HEIGHT);
    const hdr_shader = c.LoadShader(null, r_HDR_FS_PATH);
    // const hdr_tex_id = c.rlLoadTexture(null, w_WIDTH, w_HEIGHT, c.RL_PIXELFORMAT_UNCOMPRESSED_R16G16B16A16, 1);
    // const hdr_tex = c.Texture2D {
    //     .id = hdr_tex_id,
    //     .width = w_WIDTH,
    //     .height = w_HEIGHT,
    //     .format = c.RL_PIXELFORMAT_UNCOMPRESSED_R16G16B16A16,
    //     .mipmaps = 1,
    // };
    // const hdr_render_tex = c.RenderTexture2D {
    //     .id = hdr_tex_id,
    //     .depth = 1,
    //     .texture = hdr_tex,
    // };
    // defer c.rlUnloadTexture(hdr_tex_id);


    // const p = Particle { .pos = .{0, 0} };
    randomize_config(random, std.heap.c_allocator);
    generate_particle(random);
    var camera = c.Camera2D {
        .offset = .{ .x = w_WIDTH/2, .y = w_HEIGHT/2 },
        .target = .{ .x = w_WIDTH/2, .y = w_HEIGHT/2 },
        .rotation = 0,
        .zoom = 1, 
    };
    const simulation_frame_rate = 60.0;
    c.SetTargetFPS(simulation_frame_rate);
    const simulation_dt = 1.0/simulation_frame_rate;
    const a = std.heap.c_allocator;
    var arena = std.heap.ArenaAllocator.init(a);
    const blend_mode: c.BlendMode = c.BLEND_ADDITIVE;
    var measurement = FixedRingBuffer(i64, 5).init_with_value(5, 0);

    
    while (!c.WindowShouldClose()) {
        const dt = c.GetFrameTime();
        _ = arena.reset(.retain_capacity);
        if (c.IsKeyPressed(c.KEY_R)) {
            randomize_config(random, std.heap.c_allocator);
            generate_particle(random);
        }
        if (c.IsKeyPressed(c.KEY_SPACE)) {
            const timestamp = std.time.microTimestamp();
            const screenshot_name = 
                std.fmt.allocPrintSentinel(arena.allocator(), "screen_shot_{}.png", .{ timestamp }, 0) catch unreachable;
            c.TakeScreenshot(screenshot_name);
        }

        const camera_move_spd = 500;
        const camera_zoom_spd = 5;
        camera.zoom += c.GetMouseWheelMove() * camera_zoom_spd * dt;
        camera.zoom = @max(1, camera.zoom);
        if (c.IsKeyDown(c.KEY_W)) {
            camera.target.y += camera_move_spd * dt;
        }
        if (c.IsKeyDown(c.KEY_S)) {
            camera.target.y -= camera_move_spd * dt;
        }
        if (c.IsKeyDown(c.KEY_D)) {
            camera.target.x += camera_move_spd * dt;
        }
        if (c.IsKeyDown(c.KEY_A)) {
            camera.target.x -= camera_move_spd * dt;
        }
             
        const update_start = std.time.milliTimestamp();
        update_game(simulation_dt);
        const update_end = std.time.milliTimestamp();
        measurement.push(update_end - update_start);


        c.BeginTextureMode(HDRBuffer);
        {
            c.ClearBackground(c.BLACK);
            c.BeginShaderMode(shader);
            c.BeginMode2D(camera);
            c.BeginBlendMode(@intCast(blend_mode));

            for (particles) |p| {
                DrawParticle(p);
            }
            c.EndBlendMode();
            c.EndMode2D();
            c.EndShaderMode();
        }
        c.EndTextureMode();

        c.BeginDrawing();
        {

            c.BeginShaderMode(hdr_shader);
            c.DrawTexture(HDRBuffer.texture, 0, 0, c.WHITE);
            c.EndShaderMode();

            const measurement_txt = 
                std.fmt.allocPrintSentinel(
                    arena.allocator(),
                    "update time: {} ms",
                    .{ measurement.avg() }, 0) catch unreachable;
            const frame_time_txt = 
                std.fmt.allocPrintSentinel(
                    arena.allocator(),
                    "frame time: {} ms",
                    .{ dt * 1000 }, 0) catch unreachable;

            c.DrawText(measurement_txt, 50, 50, 20, c.WHITE);
            c.DrawText(frame_time_txt, 50, 75, 20, c.WHITE);
            c.DrawFPS(50, 100);
        }
        c.EndDrawing();
    }
}
