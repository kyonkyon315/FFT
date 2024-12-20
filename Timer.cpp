#include "Timer.h"


Timer::Timer()
	:
	running(true), m_stop()
{
	m_start = std::chrono::system_clock::now();
}

void Timer::start() {
	running = true;
	m_start = std::chrono::system_clock::now();
}

void Timer::stop() {
	m_stop = std::chrono::system_clock::now();
	if (!running) {
		throw std::runtime_error(
			" In void Timer::stop() \n"
			" timer is not started.\n");
	}
	running = false;
}

double Timer::getTime() {
	if (running) {
		throw std::runtime_error(
			" In double Timer::getTime() \n"
			" timer is running.\n");
	}
	auto duration=std::chrono::duration_cast<std::chrono::microseconds>(m_stop - m_start);
	return (double)duration.count();
}

std::ostream& operator<<(std::ostream& ostm, const Timer& timer) {
	if (timer.running) {
		throw std::runtime_error(
			" In std::ostream& operator<<(std::ostream& ostm, const Timer& timer) \n"
			" timer is running.\n");
	}
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(timer.m_stop - timer.m_start);
	ostm << "time " << std::setw(10) << duration.count() << "[microseconds]";
	return ostm;
}