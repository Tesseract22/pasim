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

const g_PARTICLE_TEXTURE_SIZE = 20;
const g_PARTICLE_RADIUS: f32 = 0.12;


// Assume already in shader mode
fn DrawParticle(p: Particle) void {
    const hsv: c.Vector3 = switch (p.kind) {
        0 => .{ .x = 5, .y = 0.8, .z = 1.0 },
        1 => .{ .x = 200, .y = 0.8, .z = 0.75 },
        2 => .{ .x = 100, .y = 0.8, .z = 0.75 },
        else => unreachable,
    };
    const rgb = c.ColorFromHSV(hsv.x, hsv.y, hsv.z);
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

const particle_kind = 3;
const force_config_matrix = [particle_kind][particle_kind]ForceConfig {
    .{
        .{ .radius = 0.07, .strength = 0.01 },
        .{ .radius = 0.09, .strength = -0.04 },
        .{ .radius = 0.02, .strength = -0.01 },
    },
    .{
        .{ .radius = 0.07, .strength = 0.02 },
        .{ .radius = 0.10, .strength = -0.03 },
        .{ .radius = 0.10, .strength = -0.04 },
    },
    .{
        .{ .radius = 0.05, .strength = 0.02 },
        .{ .radius = 0.05, .strength = -0.03 },
        .{ .radius = 0.05, .strength = -0.01 },
    }

};

fn linear_force(cfg: ForceConfig, d: f32) f32 {
    const f = cfg.strength * @max(0, (cfg.radius-@abs(d)) / cfg.radius);
    return f;
}

fn updateGame(dt: f32) void {
    const collision_cfg = ForceConfig {
        .radius = 0.01,
        .strength = 0.05,
    };
    for (0..particles.len) |i| {
        const a = &particles[i];
        for (i+1..particles.len) |j| {
            const b = &particles[j];
            const l = a.pos - b.pos;
            const d = dist(a.pos, b.pos);
            const unit_l = if (d == 0) splatv2(0) else l / splatv2(d);
            // if (d > collision_max_dist) continue;
            const collision_force = linear_force(collision_cfg, d) * a.mass * b.mass;

            const interact_cfg_ab = force_config_matrix[a.kind][b.kind];
            const interact_cfg_ba = force_config_matrix[b.kind][a.kind];
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
        a.spd -= a.spd * splatv2(0.2 * dt);
    }
}

pub fn main() !void {
    c.InitWindow(w_WIDTH, w_HEIGHT, "Pasim"); 
    particle_default_texture = c.Texture2D {
        .id = c.rlGetTextureIdDefault(),
        .height = 1,
        .width = 1,
        .mipmaps = 1,
        .format = c.PIXELFORMAT_UNCOMPRESSED_R8G8B8A8,
    };

    var rand = std.Random.Xoroshiro128.init(0);
    const random = rand.random();



    const shader = c.LoadShader(null, r_FS_PATH);
    const shader_radius_uniform = c.GetShaderLocation(shader, "radius");
    c.SetShaderValue(shader, shader_radius_uniform, @ptrCast(&g_PARTICLE_RADIUS), c.SHADER_UNIFORM_FLOAT);

    // const p = Particle { .pos = .{0, 0} };
    const init_vel_mul = 0.1;
    for (&particles) |*p| {
        p.pos = .{ (random.float(f32)-0.5)*2*w_RATIO, (random.float(f32)-0.5)*2 };
        p.spd = .{ (random.float(f32)-0.5)*init_vel_mul, (random.float(f32)-0.5)*init_vel_mul };
        p.mass = (random.float(f32) + 0.5) * 2;
        p.kind = random.int(u8) % particle_kind;
    }
    const camera = c.Camera2D {
        .offset = .{ .x = w_WIDTH/2, .y = w_HEIGHT/2 },
        .target = .{ .x = w_WIDTH/2, .y = w_HEIGHT/2 },
        .rotation = 0,
        .zoom = 1, 
    };
    c.SetTargetFPS(60);
    while (!c.WindowShouldClose()) {
        const dt = c.GetFrameTime();
        updateGame(dt);

        c.BeginDrawing(); 
        c.ClearBackground(c.BLACK);


        c.BeginShaderMode(shader);
        c.BeginMode2D(camera);

        // DrawParticle(p);
        for (particles) |p| {
            DrawParticle(p);
        }

        c.EndMode2D();
        c.EndShaderMode();

        c.EndDrawing();
    }
}
