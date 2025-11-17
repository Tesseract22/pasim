const std = @import("std");

const Perf = @This();

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

const Measurement = FixedRingBuffer(struct {i64, i64}, 5);
var measurements: std.StringHashMap(Measurement) = undefined;

pub fn init(a: std.mem.Allocator) void {
    measurements = .init(a);
}

pub fn measure_start(name: []const u8) void {
    const gop = measurements.getOrPut(name) catch unreachable;
    if (!gop.found_existing) {
        gop.value_ptr.* = .init_with_value(5, .{0,0});
    }
    gop.value_ptr.push(.{std.time.milliTimestamp(), undefined});
}

pub fn measure_end(name: []const u8) void {
    const measurement = measurements.getPtr(name) orelse std.process.fatal("no measurement named {s}", .{ name });
    measurement.last().*[1] = std.time.milliTimestamp();
    // std.log.err("time: {}, avg: {}", .{measurement.last().*, measurement.avg() });
}

pub fn measure_avg(comptime name: []const u8) i64 {
    return measure_avg_safe(name) orelse std.process.fatal("no measurement named {s}", .{ name });
}

pub fn measure_avg_safe(comptime name: []const u8) ?i64 {
    const measurement = measurements.getPtr(name) orelse return null;
    var avg: i64 = 0;
    for (0..5) |i| {
        avg += measurement.data[i][1] - measurement.data[i][0];
    }
    return @divFloor(avg, 5);
}
