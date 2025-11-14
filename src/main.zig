const std = @import("std");

const c = @cImport({
    @cInclude("raylib.h");
    @cInclude("rlgl.h");
});

const w_HEIGHT = 1080;
const w_WIDTH = 1920;
const w_VIEWPORT = w_HEIGHT;
const w_RATIO: comptime_float = 
    @as(comptime_float, @floatFromInt(w_WIDTH)) / 
    @as(comptime_float, @floatFromInt(w_HEIGHT));

const r_VS_PATH = "resources/base.vs";
const r_FS_PATH = "resources/base.fs";
const r_HDR_FS_PATH = "resources/hdr.fs";


 // normalized world coordinate
 // (0, 1) is the top of screen
const Vec2 = @Vector(2, f32);
const WorldCoord =  @Vector(2, f32);
const ScreenCoord =  @Vector(2, f32);
const Particle = struct {
    pos: WorldCoord,
    spd: @Vector(2, f32),
    mass: f32 = 1,
    kind: u8,
};

var particles: [60000]Particle = undefined;

fn splatv2(f: f32) @Vector(2, f32) {
    return @splat(f);
}

fn coord2rlcoord(p: WorldCoord) c.Vector2 {
    return .{ .x = p[0], .y = p[1] };
}

fn wdist2s(d: f32) f32 {
    return d/2.0 * w_VIEWPORT;
}

fn sdist2w(d: f32) f32 {
    return d/w_VIEWPORT * 2;
}

fn w2s(p: WorldCoord) ScreenCoord {
    const v = (p + splatv2(1)) / splatv2(2) * splatv2(w_VIEWPORT); // this is from 0-1
    return .{ v[0]+(w_WIDTH-w_HEIGHT)/2, w_HEIGHT-v[1] };
}

fn s2w(p: ScreenCoord) WorldCoord {
    const v = ScreenCoord { p[0]-(w_WIDTH-w_HEIGHT)/2, p[1] };
    return v / splatv2(w_VIEWPORT) * splatv2(2) - splatv2(1);
}

fn flip_y(p: @Vector(2, f32)) @Vector(2, f32) {
    return .{ p[0], -p[1] };
}

fn size2pixels(s: f32) f32 {
    return s / 2.0 * w_VIEWPORT;
}

// fn s2w_camera(p: ScreenCoord, camera: c.Camera2D) WorldCoord {
//     const camera_pos = ScreenCoord { camera.target.x, camera.target.y };
//     const center = ScreenCoord { w_WIDTH/2, w_HEIGHT/2 };
//     const d = (p - center) / splatv2(w_VIEWPORT/(2*camera.zoom));
//     return flip_y(d) + s2w(camera_pos);
// }

fn dist2(a: WorldCoord, b: WorldCoord) f32 {
    const d = a - b;
    return @reduce(.Add, d*d);
}

fn dist(a: WorldCoord, b: WorldCoord) f32 {
    return @sqrt(dist2(a, b));
}

var particle_default_texture: c.Texture2D = undefined;

const g_PARTICLE_RADIUS: f32 = 0.0025;
const g_PARTICLE_TEXTURE_SIZE: f32 = 100;
const g_PARTICLE_TEXTURE_RADIUS: f32 = g_PARTICLE_RADIUS / (g_PARTICLE_TEXTURE_SIZE / w_WIDTH);


// Assume already in shader mode
fn DrawParticle(p: Particle) void {
    const rgb = particle_colors.items[p.kind];
    const sreen_coord = w2s(.{p.pos[0], -p.pos[1]});
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

var particle_kind: u32 = 0;
var particle_force_configs = std.ArrayList(ForceConfig).empty;
var particle_colors = std.ArrayList(c.Color).empty;
const collision_cfg = ForceConfig {
    .radius = 0.025,
    .strength = 0.6,
};
const particle_drag = 100;
const r_force_radius_max = 0.099;
const r_force_radius_min = 0.05;
const r_force_strength_max = 0.20;
const r_force_strength_min = 0.10;

const mouse_force = ForceConfig {
    .radius = 0.1,
    .strength = 5,
};

const p_GRID_SIZE = 0.1;
const p_GRID_SPACING = wdist2s(p_GRID_SIZE);
const p_GRID_H_SLICES: usize = @intFromFloat(@ceil(@as(comptime_float, w_HEIGHT)/p_GRID_SPACING)); 
const p_GRID_V_SLICES: usize = @intFromFloat(@ceil(@as(comptime_float, w_WIDTH)/p_GRID_SPACING)); 
const p_GRID_CELL = p_GRID_H_SLICES * p_GRID_V_SLICES;
var grid_bins: [particles.len]u32 = undefined;
// the i'th bin contains particles from grid_bins[grid_bins_range[i] : grid_bins_range[i+1]]
var grid_bins_range: [p_GRID_CELL+1]u32 = undefined;

comptime {
    const assert = std.debug.assert;
    assert(p_GRID_SIZE > r_force_radius_max);
    assert(r_force_radius_min > collision_cfg.radius);
    assert(r_force_strength_max < collision_cfg.strength);
}

fn pos_in_bin(wpos: WorldCoord) u32 {
    const spos = w2s(wpos);
    const grid_pos: @Vector(2, u32) = @intFromFloat(spos / splatv2(p_GRID_SPACING));
    const clamped = @Vector(2, u32) { std.math.clamp(grid_pos[0], 0, @as(u32, p_GRID_V_SLICES)-1), std.math.clamp(grid_pos[1], 0, @as(u32, p_GRID_H_SLICES)-1) };
    const idx = get_grid_index(clamped[0], clamped[1]);
    return idx;
} 

fn compute_bin() void {
    var bin_sizes: [p_GRID_CELL]u32 = undefined;
    var bin_sizes_prefix: [p_GRID_CELL+1]u32 = undefined;

    @memset(&bin_sizes, 0);
    for (particles) |p| {
        const grid_index = pos_in_bin(p.pos);
        // std.log.debug("p: {}, grid: {}", .{p.pos, grid_index});
        bin_sizes[grid_index] += 1;
    }


    var sum: u32 = 0;
    for (bin_sizes, 0..) |size, i| {
        bin_sizes_prefix[i] = sum;
        sum += size; 
    }
    bin_sizes_prefix[bin_sizes_prefix.len-1] = sum;

    @memcpy(&grid_bins_range, &bin_sizes_prefix);

    for (particles, 0..) |p, i| {
        const grid_index = pos_in_bin(p.pos);
        const bin_offset = &bin_sizes_prefix[grid_index];
        grid_bins[bin_offset.*] = @intCast(i);
        bin_offset.* += 1;
    }
}

var global_random: std.Random = undefined;


fn float_range(random: std.Random, min: f32, max: f32) f32 {
    return random.float(f32) * (max - min) + min;
}

fn random_sign(random: std.Random) f32 {
    return if (random.boolean()) 1 else -1;
}

fn get_random_unit_sphere(random: std.Random) Vec2 {
    const angle = 2*std.math.pi * random.float(f32);
    return .{
        @sin(angle),
        @cos(angle),
    };
}

fn randomize_config(a: std.mem.Allocator) void {
    const random = global_random;
    particle_force_configs.clearRetainingCapacity();
    particle_colors.clearRetainingCapacity();
    // particle_kind = random.int(u8) % 5 + 2;
    particle_kind = 10;
    for (0..particle_kind) |i| {
        for (0..particle_kind) |j| {
            _ = i;
            _ = j;
            particle_force_configs.append(a, 
                .{
                    .radius = float_range(random, r_force_radius_min, r_force_radius_max),
                    .strength = random_sign(random) * float_range(random, r_force_strength_min, r_force_strength_max),
                })
            catch unreachable;
        }
        const hsv =  c.Vector3 { .x = random.float(f32) * 360, .y = float_range(random, 0.5, 1), .z = float_range(random, 0.5, 1) };
        const rgb = c.ColorFromHSV(hsv.x, hsv.y, hsv.z);
        particle_colors.append(a, rgb) catch unreachable;
    }
}

fn generate_particle() void {
    // const init_vel_mul = 0.1;
    const random = global_random;
    const init_vel_mul = 0;
    for (&particles, 0..) |*p, i| {
        p.pos = .{ (random.float(f32)-0.5)*2*w_RATIO, (random.float(f32)-0.5)*2 };
        p.spd = .{ (random.float(f32)-0.5)*init_vel_mul, (random.float(f32)-0.5)*init_vel_mul };
        // p.mass = (random.float(f32) + 0.5) * 2;
        p.mass = 1;
        p.kind = @intCast(i % particle_kind);
    }
}

fn linear_force(cfg: ForceConfig, d: f32) f32 {
    const radius = cfg.radius;
    const strength = cfg.strength;
    const f = strength * @max(0, (radius-@abs(d)) / radius);
    return f;
}

fn compute_interaction(a: *Particle, b: *Particle, dt: f32) void {
    const l = a.pos - b.pos;
    const d = dist(a.pos, b.pos); // since we uses linear_force, zero distance does not cause infinite force
    const unit_l = if (d == 0) .{1,0} else l / splatv2(d);
    // if (d > collision_max_dist) continue;
    const collision_force = linear_force(collision_cfg, d) * a.mass * b.mass;

    const interact_cfg_ab = particle_force_configs.items[a.kind * particle_kind + b.kind];
    //const interact_cfg_ba = particle_force_configs.items[b.kind * particle_kind + a.kind];
    const interact_force_ab = linear_force(interact_cfg_ab, d);
    //const interact_force_ba = linear_force(interact_cfg_ba, d);
    // c.DrawLineEx(coord2rlcoord(w2s(a.pos)), coord2rlcoord(w2s(b.pos)), 1, c.ORANGE);

    a.spd += splatv2(dt*interact_force_ab) * unit_l;
    // b.spd -= splatv2(dt*interact_force_ba) * unit_l;

    a.spd += splatv2(dt*collision_force) * unit_l;
    //b.spd -= splatv2(dt*collision_force/2) * unit_l;

    if (d == 0) {
        const random_unit = get_random_unit_sphere(global_random);
        a.pos += random_unit * splatv2(1e-3);
        b.pos -= random_unit * splatv2(1e-3);
    }

}

fn compute_interaction_in_range(particle: u32, bin_start: u32, bin_end: u32, dt: f32) void {
    const a = &particles[particle];
    for (bin_start..bin_end) |j| {
        const other = grid_bins[j]; 
        if (other == particle) continue;
        const b = &particles[other];
        compute_interaction(a, b, dt); 
    }
}

fn compute_interaction_with_bin(particle: u32, grid_index: u32, dt: f32) void {
    const bin_start = grid_bins_range[grid_index];
    const bin_end = grid_bins_range[grid_index+1];
    compute_interaction_in_range(particle, bin_start, bin_end, dt);
}


fn compute_interaction_in_bin(grid_index: u32, dt: f32) void {
    const grid_x: u32 = @intCast(grid_index % p_GRID_V_SLICES);
    const grid_y: u32 = @intCast(grid_index / p_GRID_V_SLICES);
    const bin_start = grid_bins_range[grid_index];
    const bin_end = grid_bins_range[grid_index+1];
    for (bin_start..bin_end) |i| {
        const particle: u32 = grid_bins[i];
        compute_interaction_with_bin(particle, @intCast(grid_index), dt);
        if (grid_x > 0) compute_interaction_with_bin(particle,
            get_grid_index(grid_x-1, grid_y), dt);
        if (grid_x < p_GRID_V_SLICES-1) compute_interaction_with_bin(particle,
            get_grid_index(grid_x+1, grid_y), dt);
        if (grid_y > 0) compute_interaction_with_bin(particle,
            get_grid_index(grid_x, grid_y-1), dt);
        if (grid_y < p_GRID_H_SLICES-1) compute_interaction_with_bin(particle,
            get_grid_index(grid_x, grid_y+1), dt);

        if (grid_x > 0 and
            grid_y > 0) 
            compute_interaction_with_bin(particle,
                get_grid_index(grid_x-1, grid_y-1), dt);
        if (grid_x > 0 and
            grid_y < p_GRID_H_SLICES-1)
            compute_interaction_with_bin(particle,
                get_grid_index(grid_x-1, grid_y+1), dt);
        if (grid_x < p_GRID_V_SLICES-1 and
            grid_y > 0)
            compute_interaction_with_bin(particle,
                get_grid_index(grid_x+1, grid_y-1), dt);
        if (grid_x < p_GRID_V_SLICES-1 and
            grid_y < p_GRID_H_SLICES-1)
            compute_interaction_with_bin(particle,
                get_grid_index(grid_x+1, grid_y+1), dt);
    }
}

fn compute_interaction_in_bin_parallel(grid_index_a: *std.atomic.Value(u32), dt: f32) void {
    while (true) {
        const grid_index = grid_index_a.fetchAdd(1, .monotonic);
        if (grid_index >= p_GRID_CELL) break;
        compute_interaction_in_bin(grid_index, dt);
    }
}

fn get_grid_index(x: u32, y: u32) u32 {
    return x + y * @as(u32, p_GRID_V_SLICES);
}

var collision_method_basic = false;


fn update_game(dt: f32, mouse_pos: WorldCoord, mouse_action: ?enum { attract, repel }) void {
    if (collision_method_basic) {
        for (0..particles.len) |i| {
            const a = &particles[i];
            for (i+1..particles.len) |j| {
                const b = &particles[j];
                compute_interaction(a, b, dt); 
            }
        }
    } else {
        compute_bin();
        var threads: [10]std.Thread = undefined;
        var grid_index = std.atomic.Value(u32).init(0);
        for (&threads) |*t| {
            t.* = std.Thread.spawn(.{}, compute_interaction_in_bin_parallel, .{ &grid_index, dt }) catch unreachable;
        }
        for (&threads) |*t| {
            t.join();
        }
        // for (0..p_GRID_CELL) |grid_index| {
        //     compute_interaction_in_bin(@intCast(grid_index), dt);
        // }
    }

    for (0..particles.len) |i| {
        const a = &particles[i];
        a.spd *= splatv2(@exp(-particle_drag*dt*dist(a.spd, .{0,0})));
        // mouse force
        const d = dist(a.pos, mouse_pos);
        const l = a.pos - mouse_pos;
        const unit_l = if (d == 0) splatv2(0) else l / splatv2(d);
        if (mouse_action) |action| {
            const f = linear_force(mouse_force, d);
            switch (action) {
                .attract =>
                    a.spd += splatv2(dt*f/a.mass) * unit_l,
                .repel =>
                    a.spd -= splatv2(dt*f/a.mass) * unit_l,
            } 
        }

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

fn exp_smooth(current: f32, target: f32, delta: f32) f32 {
    return current + (target - current) * (1 - @exp(-delta));
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
    global_random = rand.random();



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
    randomize_config(std.heap.c_allocator);
    generate_particle();
    var camera = c.Camera2D {
        .zoom = 1, 
    };
    var camera_target = c.Camera2D {
        .zoom = 1,
    };
    const simulation_frame_rate = 30.0;
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
            randomize_config(std.heap.c_allocator);
            generate_particle();
        }
        if (c.IsKeyPressed(c.KEY_SPACE)) {
            const timestamp = std.time.microTimestamp();
            const screenshot_name = 
                std.fmt.allocPrintSentinel(arena.allocator(), "screen_shot_{}.png", .{ timestamp }, 0) catch unreachable;
            c.TakeScreenshot(screenshot_name);
        }

        // const y_shift = camera.target.y - w_HEIGHT/2;
        const mouse_abs_pos = c.Vector2 { .x = c.GetMousePosition().x, .y = w_HEIGHT-c.GetMousePosition().y };
        const mouse_spos = c.GetScreenToWorld2D(mouse_abs_pos, camera);
        const mouse_wpos = s2w(.{ mouse_spos.x, mouse_spos.y });
        const mouse_grid_index = pos_in_bin(mouse_wpos);

        // 
        const camera_move_spd = 500;
        const camera_zoom_spd = 0.3;
        const wheel = c.GetMouseWheelMove();
        if (wheel != 0)
        {
            camera_target.offset = mouse_abs_pos;
            camera.offset = mouse_abs_pos;
            camera_target.target = mouse_spos;
            camera.target = mouse_spos;

            // Zoom increment
            // Uses log scaling to provide consistent zoom speed
            camera_target.zoom = std.math.clamp(@exp2(@log2(camera.zoom)+wheel*camera_zoom_spd), 0.125, 64.0);
        }
        if (c.IsKeyDown(c.KEY_W)) {
            camera_target.target.y += camera_move_spd * dt;
        }
        if (c.IsKeyDown(c.KEY_S)) {
            camera_target.target.y -= camera_move_spd * dt;
        }
        if (c.IsKeyDown(c.KEY_D)) {
            camera_target.target.x += camera_move_spd * dt;
        }
        if (c.IsKeyDown(c.KEY_A)) {
            camera_target.target.x -= camera_move_spd * dt;
        }
        if (c.IsKeyPressed(c.KEY_Z)) {
            collision_method_basic = !collision_method_basic;
        }
        const camera_smooth_spd = 10;
        camera.zoom = exp_smooth(camera.zoom, camera_target.zoom, dt * camera_smooth_spd);
        camera.target.x = exp_smooth(camera.target.x, camera_target.target.x, dt * camera_smooth_spd);
        camera.target.y = exp_smooth(camera.target.y, camera_target.target.y, dt * camera_smooth_spd);
        camera.offset.x = exp_smooth(camera.offset.x, camera_target.offset.x, dt * camera_smooth_spd);
        camera.offset.y = exp_smooth(camera.offset.y, camera_target.offset.y, dt * camera_smooth_spd);

        

        const update_start = std.time.milliTimestamp();
        update_game(simulation_dt, mouse_wpos, if (c.IsMouseButtonDown(c.MOUSE_BUTTON_LEFT)) .attract else if (c.IsMouseButtonDown(c.MOUSE_BUTTON_RIGHT)) .repel else null);
        const update_end = std.time.milliTimestamp();
        measurement.push(update_end - update_start);

        const mouse_spos_test = w2s(.{ mouse_wpos[0], -mouse_wpos[1] }) ;
        // std.log.err("zoom: {}, mouse: spos: {} wpos: {}", .{ camera.zoom, mouse_spos, mouse_wpos});
        

        c.BeginTextureMode(HDRBuffer);
        {
            c.ClearBackground(c.BLACK);
            c.BeginMode2D(camera);

            c.DrawCircle(@intFromFloat(mouse_spos_test[0]), @intFromFloat(mouse_spos_test[1]), size2pixels(mouse_force.radius),
                .{.r = 0xff, .g = 0xff, .b = 0xff, .a = 0x3f });

            c.BeginBlendMode(@intCast(blend_mode));
            c.BeginShaderMode(shader);
            for (particles) |p| {
                DrawParticle(p);
            }
            c.EndShaderMode();
            c.EndBlendMode();
        
            if (!collision_method_basic and false) {
                for (0..p_GRID_H_SLICES+1) |y| {
                    const yf: f32 = @floatFromInt(y);
                    c.DrawLineEx(.{.x=0, .y=yf*p_GRID_SPACING }, .{.x=w_WIDTH, .y=yf*p_GRID_SPACING}, 1, c.WHITE);
                }
                for (0..p_GRID_V_SLICES+1) |x| {
                    const xf: f32 = @floatFromInt(x);
                    c.DrawLineEx(.{.x=xf*p_GRID_SPACING, .y=0}, .{.x=xf*p_GRID_SPACING, .y=w_HEIGHT}, 1, c.WHITE);
                }
                const grid_x = mouse_grid_index % p_GRID_V_SLICES;
                const grid_y = mouse_grid_index / p_GRID_V_SLICES;
                const grid_index = get_grid_index(@intCast(grid_x), @intCast(grid_y));
                const grid_start = grid_bins_range[grid_index];
                const grid_end = grid_bins_range[grid_index+1];

                std.log.debug("mouse at pos: {}, grid: {},{} [{}], {} particles", .{ mouse_wpos, grid_x, grid_y, grid_index, grid_end-grid_start });
                //for (grid_start..grid_end) |i| {
                //    const pi = grid_bins[i];
                //    const p = particles[pi];
                //    std.log.debug("p: {}, {}", .{ pi, p.pos });
                //}
            }

            c.EndMode2D();
        }
        c.EndTextureMode();

        c.BeginDrawing();
        {

            c.ClearBackground(c.BLACK);
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
